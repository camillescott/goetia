/* compactor.hh -- streaming dBG compactor
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef COMPACTOR_HH
#define COMPACTOR_HH

#include <assert.h>

#include "boink/assembly.hh"
#include "boink/hashing.hh"
#include "boink/dbg.hh"
#include "boink/cdbg/cdbg.hh"
#include "boink/minimizers.hh"
#include "boink/event_types.hh"
#include "boink/reporting/report_types.hh"


namespace boink {
namespace cdbg {

using namespace boink::event_types;
using namespace boink::reporting::report_types;

# ifdef DEBUG_CPTR
#   define pdebug(x) do { std::cerr << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif

#define complement(ch) ((ch) == 'A' ? 'T' : \
                        (ch) == 'T' ? 'A' : \
                        (ch) == 'C' ? 'G' : 'C')


/* Represents a segment of new k-mers from a sequence, relative to the current
 * state of the cDBG. Can either be: a null segment representing a portion
 * of the sequence already in the graph; a decision k-mer; or a new unitig or
 * portion of a unitig.
 */

struct compact_segment {
    hash_t left_anchor;
    hash_t right_anchor;
    hash_t left_flank;
    hash_t right_flank;

    // whether this segment represents a decision k-mer
    bool is_decision_kmer;
    // start position of the segment within the originating sequence
    size_t start_pos;
    // length of the segment sequence (from beginning of first k-mer to
    // end of last k-mer)
    size_t length;
    // tags associated with this segment
    HashVector tags;

    // the default constructor creates a null segment
    compact_segment()
        : left_anchor(0),
          right_anchor(0),
          left_flank(0),
          right_flank(0),
          is_decision_kmer(false),
          start_pos(0),
          length(0) {}

    compact_segment(hash_t left_anchor,
                    hash_t right_anchor,
                    bool is_decision_kmer,
                    size_t start_pos,
                    size_t length)
        : left_anchor(left_anchor),
          right_anchor(right_anchor),
          is_decision_kmer(is_decision_kmer),
          start_pos(start_pos),
          length(length) {}

    // return if the segment is default constructed / null
    // used as a delimiter token between connected runs of segments
    const bool is_null() const {
        return length == 0 && (left_anchor == right_anchor) && !is_decision_kmer;
    }
};


std::ostream& operator<<(std::ostream& os, const compact_segment& segment)
{
    os << "<compact_segment left flank=" << segment.left_flank
       << " left_anchor=" << segment.left_anchor
       << " right_anchor=" << segment.right_anchor
       << " right_flank=" << segment.right_flank
       << " start=" << segment.start_pos
       << " length=" << segment.length
       << ">";
    return os;
}


template <class GraphType>
class StreamingCompactor : public CompactorMixin<GraphType>,
                           public EventNotifier {

protected:

    uint64_t _minimizer_window_size;
    bool _cdbg_external;

public:

    using ShifterType = typename GraphType::shifter_type;
    using CompactorType = CompactorMixin<GraphType>;
    using MinimizerType = WKMinimizer<ShifterType>;
    using cDBGType = cDBG<GraphType>;

    using CompactorType::filter_nodes;
    using CompactorType::find_left_kmers;
    using CompactorType::find_right_kmers;
    using CompactorType::gather_left;
    using CompactorType::gather_right;
    using CompactorType::compactify_left;
    using CompactorType::compactify_right;
    using CompactorType::get_decision_neighbors;

    typedef ShifterType shifter_type;
    typedef GraphType graph_type;
    typedef MinimizerType minimizer_type;
    typedef cDBGType cdbg_type;

    shared_ptr<GraphType> dbg;
    shared_ptr<cDBG<GraphType>> cdbg;

    StreamingCompactor(shared_ptr<GraphType> dbg,
                       uint64_t minimizer_window_size=8,
                       shared_ptr<cDBG<GraphType>> cdbg=nullptr)
        : CompactorMixin<GraphType>(dbg),
          EventNotifier(),
          _minimizer_window_size(minimizer_window_size),
          dbg(dbg)
    {
        if (cdbg == nullptr) {
            this->cdbg = make_shared<cDBG<GraphType>>(dbg);
            _cdbg_external = false;
        } else {
            this->cdbg = cdbg;
            _cdbg_external = true;
        }
    }

    ~StreamingCompactor() {
        //cdbg->wait_on_processing(0);

        // make sure nothing else has a lock on the cdbg
        auto lock = cdbg->lock_nodes();
    }
    
    StreamingCompactorReport get_report() {
        StreamingCompactorReport report;
        report.n_full            = cdbg->metrics->n_full.Value();
        report.n_tips            = cdbg->metrics->n_tips.Value();
        report.n_islands         = cdbg->metrics->n_islands.Value();
        report.n_trivial         = cdbg->metrics->n_trivial.Value();
        report.n_circular        = cdbg->metrics->n_circular.Value();
        report.n_loops           = cdbg->metrics->n_loops.Value();
        report.n_dnodes          = cdbg->metrics->n_dnodes.Value();
        report.n_unodes          = cdbg->metrics->n_unodes.Value();
        report.n_tags            = cdbg->n_tags();
        report.n_updates         = cdbg->n_updates();
        report.n_splits          = cdbg->metrics->n_splits.Value();
        report.n_merges          = cdbg->metrics->n_merges.Value();
        report.n_extends         = cdbg->metrics->n_extends.Value();
        report.n_clips           = cdbg->metrics->n_clips.Value();
        report.n_deletes         = cdbg->metrics->n_deletes.Value();
        report.n_circular_merges = cdbg->metrics->n_circular_merges.Value();
        report.n_unique          = dbg->n_unique();
        report.estimated_fp      = dbg->estimated_fp();

        return report;
    }

    void update_sequence(const string& sequence) {
        set<hash_t> new_kmers;
        deque<compact_segment> segments;
        set<hash_t> new_decision_kmers;
        deque<NeighborBundle> decision_neighbors;
        vector<hash_t> hashes;

        find_new_segments(sequence,
                          hashes,
                          new_kmers,
                          segments,
                          new_decision_kmers,
                          decision_neighbors);

        update_from_segments(sequence,
                             new_kmers,
                             segments,
                             new_decision_kmers,
                             decision_neighbors);

        for (auto h : hashes) {
            dbg->add(h);
        }
    }

    compact_segment init_segment(hash_t left_anchor,
                                 hash_t left_flank,
                                 size_t start_pos) {
        compact_segment segment;
        segment.start_pos = start_pos;
        segment.left_anchor = left_anchor;
        segment.left_flank = left_flank;
        segment.is_decision_kmer = false;
        pdebug("Init segment: start=" << start_pos
               << " left_flank=" << left_flank
               << " left_anchor=" << left_anchor);

        return segment;
    }

    void finish_segment(compact_segment& segment,
                        size_t end,
                        hash_t right_anchor,
                        hash_t right_flank,
                        deque<compact_segment>& segments) {
        
        segment.length = end - segment.start_pos + this->_K;
        segment.right_anchor = right_anchor;
        segment.right_flank = right_flank;
        segments.push_back(segment);
        pdebug("Finished segment: " << segment);
    }

    void finish_decision_segment(compact_segment& segment,
                                 deque<compact_segment>& segments) {
        segment.length = this->_K;
        segment.right_anchor = segment.left_anchor;
        segment.right_flank = segment.left_flank;
        segment.is_decision_kmer = true;
        segments.push_back(segment);

        pdebug("Finished decision segment: " << segment);       
    }

    void find_new_segments(const string& sequence,
                           deque<compact_segment>& result) {
        // convenience function for cython land
        vector<hash_t> hashes;
        set<hash_t> new_kmers;
        set<hash_t> new_decision_kmers;
        deque<NeighborBundle> decision_neighbors;

        find_new_segments(sequence,
                          hashes,
                          new_kmers,
                          result,
                          new_decision_kmers,
                          decision_neighbors);

    }

    void find_new_segments(const string& sequence,
                           vector<hash_t>& hashes,
                           set<hash_t>& new_kmers,
                           deque<compact_segment>& segments,
                           set<hash_t>& new_decision_kmers,
                           deque<NeighborBundle>& decision_neighbors) {

        pdebug("FIND SEGMENTS: " << sequence);

        KmerIterator<ShifterType> kmers(sequence, this->_K);
        hash_t prev_hash, cur_hash;
        size_t pos = 0;
        bool cur_new = false, prev_new = false, cur_seen = false, prev_seen = false;

        deque<compact_segment> preprocess;
#ifdef DEBUG_CPTR
        vector<count_t> counts;
#endif
        compact_segment current_segment; // start null
        vector<CompactorType> segment_shifters;
        while(!kmers.done()) {
            cur_hash = kmers.next();
            cur_new = this->dbg->get(cur_hash) == 0;
            cur_seen = new_kmers.count(cur_hash);
            hashes.push_back(cur_hash);
#ifdef DEBUG_CPTR
            counts.push_back(this->dbg->get(cur_hash));
#endif

            if(cur_new && !cur_seen) {

                if(!prev_new || prev_seen) {
                    pdebug("old -> new (pos=" << pos << ")");

                    preprocess.push_back(current_segment);
                    current_segment = init_segment(cur_hash,
                                                   prev_hash,
                                                   pos);
                    segment_shifters.emplace_back(CompactorType(dbg,
                                                                *(kmers.shifter)));
                }
                new_kmers.insert(cur_hash);
            } else if (!current_segment.is_null() &&
                        (!cur_new || cur_seen)) {
                pdebug("new -> old, was_seen=" << cur_seen);
                finish_segment(current_segment, 
                               pos - 1,
                               prev_hash,
                               cur_hash,
                               preprocess);
                current_segment = compact_segment();
            } 

            ++pos;
            prev_hash = cur_hash;
            prev_new = cur_new;
            prev_seen = cur_seen;
        }
        // make sure we close current segment if necessary
        if (cur_new && !cur_seen) {
            pdebug("sequence ended on new k-mer");
            hash_t right_flank = cur_hash;
            vector<shift_t> rneighbors = filter_nodes(kmers.shifter->gather_right(),
                                                      new_kmers);
            if (rneighbors.size() == 1) {
                right_flank = rneighbors.front().hash;
            }
            finish_segment(current_segment,
                           pos - 1, // we incr'd pos...
                           cur_hash,
                           right_flank,
                           preprocess);
        }
        preprocess.push_back(compact_segment()); // null segment

        if (preprocess.size() < 3) {
            // no new sequence
            return;
        }

        // handle edge case for left_flank of first segment if it starts at pos 0
        if (preprocess[1].start_pos == 0) {
            vector<shift_t> lneighbors = filter_nodes(segment_shifters[0].gather_left(),
                                                      new_kmers);
            if (lneighbors.size() == 1) {
                preprocess[1].left_flank = lneighbors.front().hash;
            } else {
                preprocess[1].left_flank = preprocess[1].left_anchor;
            }
        }

        // run back through segments and check for new decision kmers
        // This is all super verbose but the goal is to avoid unecessary hashing
        // and a bunch of convoluted logic up in the main segment loop
        auto shifter_iter = segment_shifters.begin();
        for (auto segment : preprocess) {
            if (segment.is_null()) {
                segments.push_back(segment);
                continue;
            }

            deque<compact_segment> decision_segments;
            size_t pos = segment.start_pos;
            size_t suffix_pos = pos + this->_K - 1;
            while (1) {
                NeighborBundle neighbors;
                if (shifter_iter->get_decision_neighbors(neighbors,
                                                         new_kmers)) {
                    decision_neighbors.push_back(neighbors);
                    new_decision_kmers.insert(shifter_iter->get());

                    auto decision_segment = init_segment(shifter_iter->get(),
                                                         shifter_iter->get(),
                                                         pos);
                    finish_decision_segment(decision_segment,
                                            decision_segments);
                }

                ++pos;
                ++suffix_pos;
                if (suffix_pos == segment.start_pos + segment.length) {
                    break;
                } else {
                    shifter_iter->shift_right(sequence[suffix_pos]);
                }
            }

            if (decision_segments.size() == 0) {
                segments.push_back(segment);
            } else if (decision_segments.size() == 1) {
                compact_segment rsegment;
                rsegment.left_flank = decision_segments.front().right_anchor;
                rsegment.left_anchor = hashes[decision_segments.front().start_pos + 1];
                rsegment.right_anchor = segment.right_anchor;
                rsegment.right_flank = segment.right_flank;
                rsegment.is_decision_kmer = false;
                rsegment.start_pos = decision_segments.front().start_pos + 1;
                rsegment.length = segment.length - (rsegment.start_pos - segment.start_pos);

                segment.right_anchor = hashes[decision_segments.front().start_pos - 1];
                segment.right_flank = decision_segments.front().left_anchor;
                segment.length = (decision_segments.front().start_pos + this->_K - 1)
                                  - segment.start_pos;

                if (segment.length >= this->_K) segments.push_back(segment);
                segments.push_back(decision_segments.front());
                if (rsegment.length >= this->_K) segments.push_back(rsegment);
            } else {
                compact_segment first;
                first.left_flank = segment.left_flank;
                first.left_anchor = segment.left_anchor;
                first.right_anchor = hashes[decision_segments.front().start_pos - 1];
                first.right_flank = hashes[decision_segments.front().start_pos];
                first.is_decision_kmer = false;
                first.start_pos = segment.start_pos;
                first.length = decision_segments.front().start_pos + this->_K - 1 - segment.start_pos;
                if (first.length >= this->_K) segments.push_back(first);
                segments.push_back(decision_segments.front());

                auto u_iter = decision_segments.begin();
                auto v_iter = decision_segments.begin() + 1;
                while (v_iter != decision_segments.end()) {
                    compact_segment new_segment;
                    new_segment.left_flank = u_iter->left_anchor;
                    new_segment.left_anchor = hashes[u_iter->start_pos+1];
                    new_segment.right_anchor = hashes[v_iter->start_pos-1];
                    new_segment.right_flank = v_iter->left_anchor;
                    new_segment.is_decision_kmer = false;
                    new_segment.start_pos = u_iter->start_pos + 1;
                    new_segment.length = (v_iter->start_pos + this->_K - 1) - new_segment.start_pos;
                    if (new_segment.length >= this->_K) segments.push_back(new_segment);
                    segments.push_back(*v_iter);
                    ++u_iter;
                    ++v_iter;
                }

                compact_segment last;
                last.right_anchor = segment.right_anchor;
                last.right_flank = segment.right_flank;
                last.left_flank = decision_segments.back().right_anchor;
                last.left_anchor = hashes[decision_segments.back().start_pos + 1];
                last.is_decision_kmer = false;
                last.start_pos = decision_segments.back().start_pos + 1;
                last.length = segment.length - (last.start_pos - segment.start_pos);
                if (last.length >= this->_K) segments.push_back(last);
            }

            ++shifter_iter;
        }

#ifdef DEBUG_CPTR
        std::ostringstream os;
        os << "k-mers: [";
        for (size_t i = 0; i < counts.size(); ++i) {
            os << std::to_string(counts[i]) << ":" << hashes[i] << ",";
        }
        os << "]";
        pdebug("FIND SEGMENTS complete: " << os.str()
                << std::endl <<  segments);
#endif
    }

   void update_from_segments(const string& sequence,
                              set<hash_t>& new_kmers,
                              deque<compact_segment>& segments,
                              set<hash_t> & new_decision_kmers,
                              deque<NeighborBundle>& decision_neighbors
                              ) {

        if (segments.size() < 3) {
            pdebug("No segments.");
            return;
        }

        pdebug(segments.size() << " segments.");

        // First find all induced decision k-kmers
        // Have to wait until all are found to proceed with induction process
        // in the cDBG though
        deque<DecisionKmer> induced;
        size_t i = 2;
        while (i < segments.size()) {
            auto u = segments.at(i-2);
            auto v = segments.at(i-1);
            auto w = segments.at(i);
            if (v.is_null()) {
                pdebug("Segment is null");
                ++i;
                continue;
            } else if (v.is_decision_kmer) {
                pdebug("Segment is decision k-mer.");
                kmer_t decision_kmer(v.left_anchor,
                                     sequence.substr(v.start_pos, this->_K));
                _build_dnode(decision_kmer);
                _find_induced_decision_nodes(decision_kmer,
                                             decision_neighbors.front(),
                                             new_kmers,
                                             induced);
                decision_neighbors.pop_front();
            } else {
                pdebug("Segment is unitig substring.");
                // v is a regular segment
                // if u is null, then we need to check for left induced d-nodes
                // or for a unitig connection
                if (u.is_null()) {

                    pdebug("Checking (u,v,w), u is null, left induce v...");
                    kmer_t root(v.left_anchor,
                                sequence.substr(v.start_pos, this->_K));
                    // TODO: future optimization. we know root is not a decision
                    // k-mer, so in this case, unless v starts at pos = 0 in the sequence,
                    // the only neighbor is the preceeding k-mer in the sequence
                    _find_induced_decision_nodes_left(root, new_kmers, induced);

                }

                // v is regular segment: if w is null, check for right induced
                // d-nodes or for a unitig connection
                if (w.is_null()) {

                    pdebug("Checking (u,v,w), w is null, right induce v...");
                    kmer_t root(v.right_anchor,
                                sequence.substr(v.start_pos + v.length - this->_K, this->_K));
                    // TODO: future optimization. we know root is not a decision
                    // k-mer, so in this case, unless v starts at pos = 0 in the sequence,
                    // the only neighbor is the preceeding k-mer in the sequence
                    _find_induced_decision_nodes_right(root, new_kmers, induced);
                }
            }
            ++i;
        }

        // Induce all the decision k-mers we found
        _induce_decision_nodes(induced, new_kmers);

        // Now, with the cDBG in a correct state, update its unitigs
        // from our new segments
        for (auto segment : segments) {
            if (!segment.is_decision_kmer && !segment.is_null()) {
                _update_unode(segment, sequence);
            }
        }
    }

    void _induce_decision_nodes(deque<DecisionKmer>& induced_decision_kmers,
                                set<hash_t>& new_kmers) {

        set<hash_t> induced_decision_kmer_hashes;

        pdebug("Perform induction on " << induced_decision_kmers.size() <<
               " new decision k-mers");
        for (auto dkmer : induced_decision_kmers) {
            _build_dnode(dkmer.first);
            induced_decision_kmer_hashes.insert(dkmer.first.hash);
        }

        set<hash_t> processed;
        size_t n_attempts = 0;
        size_t max_attempts = 4 * induced_decision_kmer_hashes.size();
        while (induced_decision_kmers.size() > 0) {
            pdebug(induced_decision_kmers.size() << " more splits to attempt...");
            n_attempts++;

            DecisionKmer d_kmer = induced_decision_kmers.front();
            induced_decision_kmers.pop_front();

            if (processed.count(d_kmer.first.hash)) {
                pdebug("Processed " << d_kmer.first << " already");
                continue;
            }

            if (_try_split_unode(d_kmer.first,
                                 d_kmer.second, 
                                 new_kmers,
                                 induced_decision_kmer_hashes,
                                 processed)) {
                pdebug("Split successful on " << d_kmer.first);
                processed.insert(d_kmer.first.hash);
            } else {
                induced_decision_kmers.push_back(d_kmer);
            }
            
            if (n_attempts > max_attempts) {
                throw BoinkException("Stuck in split attempt loop, failing.");
            }
        }
    }

    virtual bool _try_split_unode(kmer_t root,
                              NeighborBundle& neighbors,
                              set<hash_t>& new_kmers,
                              set<hash_t>& induced_decision_kmer_hashes,
                              set<hash_t>& processed) {
        pdebug("Attempt unitig split from " << root);

        UnitigNode * unode_to_split;
        
        if ((unode_to_split = cdbg->query_unode_end(root.hash)) != nullptr) {
            // special case: induced an end k-mer, just have to trim the u-node,
            // no need to create a new one

            if (unode_to_split->meta() == TRIVIAL) {
                pdebug("Induced a trivial u-node, delete it.");
                cdbg->delete_unode(unode_to_split);
                return true;
            }

            if (unode_to_split->meta() == CIRCULAR) {
                auto unitig = unode_to_split->sequence;
                cdbg->split_unode(unode_to_split->node_id,
                                  0,
                                  root.kmer,
                                  this->hash(unitig.substr(unitig.size() - this->_K)),
                                  this->hash(unitig.c_str() +1));
                return true;
            }

            hash_t new_end;
            direction_t clip_from;
            if (root.hash == unode_to_split->left_end()) {
                new_end = this->hash(unode_to_split->sequence.c_str() + 1);
                clip_from = DIR_LEFT;
            } else {
                new_end = this->hash(unode_to_split->sequence.c_str()
                                     + unode_to_split->sequence.size()
                                     - this->_K - 1);
                clip_from = DIR_RIGHT;
            }
            cdbg->clip_unode(clip_from,
                             root.hash,
                             new_end);
            return true;
        }

        vector<kmer_t> lfiltered;
        std::copy_if(neighbors.first.begin(),
                     neighbors.first.end(),
                     std::back_inserter(lfiltered),
                     [&] (kmer_t neighbor) { return
                        !new_kmers.count(neighbor.hash) &&
                        !processed.count(neighbor.hash);
                     });

        vector<kmer_t> rfiltered;
        std::copy_if(neighbors.second.begin(),
                     neighbors.second.end(),
                     std::back_inserter(rfiltered),
                     [&] (kmer_t neighbor) { return
                        !new_kmers.count(neighbor.hash) &&
                        !processed.count(neighbor.hash);
                     });
        // EDGE CASE: split k-mer has no unitig to split because of previous processed
        // nodes or because neighbor is a properly oriented unitig end

        pdebug(lfiltered.size() << " left, " << rfiltered.size() << " right");
        if (lfiltered.size()) {
            // size should always be 1 here
            pdebug("Found a valid left neighbor, search this way... ("
                   << lfiltered.size() << " in filtered set, should always be 1.)");
            auto start = lfiltered.back();
            this->set_cursor(start.kmer);
            Path path;
            hash_t end_hash;
            compactify_left(path, end_hash, processed);

            unode_to_split = cdbg->query_unode_end(end_hash);
            if (unode_to_split != nullptr) {
                size_t split_point = path.size() + 1;
                hash_t left_unode_new_right = start.hash;
                pdebug("split point is " << split_point <<
                        " new_right is " << left_unode_new_right
                       << " root was " << root.hash);
                hash_t right_unode_new_left = this->hash(unode_to_split->sequence.c_str() + 
                                                         split_point + 1);
                if (rfiltered.size()) {
                    assert(right_unode_new_left == rfiltered.back().hash);
                }

                cdbg->split_unode(unode_to_split->node_id,
                                  split_point,
                                  root.kmer,
                                  left_unode_new_right,
                                  right_unode_new_left);

                return true;
            } else {
                pdebug("Unitig to split is a loop, traversed " << this->seen.size());
                for (auto hash : this->seen) {
                    unode_to_split = cdbg->query_unode_end(hash);
                    if (unode_to_split != nullptr) break;
                }
                if (unode_to_split != nullptr) {
                    pdebug("u-node to split: " << *unode_to_split);
                    cdbg->split_unode(unode_to_split->node_id,
                                      0,
                                      root.kmer,
                                      lfiltered.back().hash,
                                      rfiltered.back().hash);
                    return true;
                }
            }
        }

        if (rfiltered.size()) {
            // size should always be 1 here
            pdebug("Found a valid right neighbor, search this way... ("
                   << rfiltered.size() << " in filtered set, should be 1.");
            auto start = rfiltered.back();
            this->set_cursor(start.kmer);
            Path path;
            hash_t end_hash;
            compactify_right(path, end_hash, processed);

            unode_to_split = cdbg->query_unode_end(end_hash);
            if (unode_to_split != nullptr) {
                size_t split_point = unode_to_split->sequence.size()
                                                     - path.size()
                                                     - 1 - this->_K;
                hash_t new_right = this->hash(unode_to_split->sequence.c_str() + 
                                              split_point - 1);
                hash_t new_left = start.hash;
                if (lfiltered.size()) {
                    assert(lfiltered.back().hash == new_right);
                }

                cdbg->split_unode(unode_to_split->node_id,
                                  split_point,
                                  root.kmer,
                                  new_right,
                                  new_left);

                return true;
            } else {
                pdebug("Unitig to split is a loop, traversed " << this->seen.size());
                for (auto hash : this->seen) {
                    unode_to_split = cdbg->query_unode_end(hash);
                    if (unode_to_split != nullptr) break;
                }
                if (unode_to_split != nullptr) {
                    pdebug("u-node to split: " << *unode_to_split);
                    cdbg->split_unode(unode_to_split->node_id,
                                      0,
                                      root.kmer,
                                      lfiltered.back().hash,
                                      rfiltered.back().hash);
                    return true;
                }
            }
        }

        // failed to find a split point. must be flanked by induced d-nodes on either
        // side before a tag or unode end: aka same unode split by three or more
        // induced d-nodes
        pdebug("Attempt to split unode failed from " << root);

        return false;
    }

    virtual void _update_unode(compact_segment& segment,
                               const string& sequence) {

        pdebug("Update Unode from segment: " << segment);

        if (segment.left_anchor == segment.right_flank
            && segment.right_flank == segment.left_anchor
            && segment.length > this->_K) {

            // TODO: What about single k-mer loop unitigs like "AAA...AAA" and so forth?
            
            pdebug("Special case: segment is a loop.");
            cdbg->build_unode(sequence.substr(segment.start_pos, segment.length),
                              segment.tags,
                              segment.left_anchor,
                              segment.left_anchor);
            return;
        }

        bool has_left_unode = cdbg->has_unode_end(segment.left_flank);
        bool has_right_unode = cdbg->has_unode_end(segment.right_flank);
        
        if (has_left_unode && !has_right_unode) {
            auto trimmed_seq = sequence.substr(segment.start_pos + this->_K - 1,
                                               segment.length - this->_K + 1);
            cdbg->extend_unode(DIR_RIGHT,
                               trimmed_seq,
                               segment.left_flank,
                               segment.right_anchor,
                               segment.tags);

        } else if (!has_left_unode && has_right_unode) {
            auto trimmed_seq = sequence.substr(segment.start_pos,
                                               segment.length - this->_K + 1);
            cdbg->extend_unode(DIR_LEFT,
                               trimmed_seq,
                               segment.right_flank,
                               segment.left_anchor,
                               segment.tags);
        } else if (has_left_unode && has_right_unode) {
            string trimmed_seq;
            pdebug("Segment is " << segment.length);
            //if (segment.length  < (this->_K * 2 - 2)) {
            //    trimmed_seq = "";
            //} else {
            //    trimmed_seq = sequence.substr(segment.start_pos + this->_K - 1,
            //                                  segment.length - (this->_K * 2 - 2));
            //}
            trimmed_seq = sequence.substr(segment.start_pos, segment.length);
            size_t span_kmers = segment.length - this->_K + 1;
            pdebug(span_kmers << " kmers in segment");

            cdbg->merge_unodes(trimmed_seq,
                               span_kmers,
                               segment.left_flank,
                               segment.right_flank,
                               segment.tags);
        } else {
            cdbg->build_unode(sequence.substr(segment.start_pos, segment.length),
                              segment.tags,
                              segment.left_anchor,
                              segment.right_anchor);
        }
    }

    uint8_t _find_induced_decision_nodes(kmer_t kmer,
                                      NeighborBundle& neighbors,
                                      set<hash_t>& new_kmers,
                                      deque<DecisionKmer>& induced) {

        return _find_induced_decision_nodes_left(kmer, neighbors, new_kmers, induced) +
               _find_induced_decision_nodes_right(kmer, neighbors, new_kmers, induced);
    }

    uint8_t _find_induced_decision_nodes_left(kmer_t kmer,
                                              set<hash_t>& new_kmers,
                                              deque<DecisionKmer>& induced) {

        pdebug("Prepare to attempt left induction on " << kmer);
        this->set_cursor(kmer.kmer);
        NeighborBundle bundle;
        bundle.first = find_left_kmers();

        if (bundle.first.size()) {
            return _find_induced_decision_nodes_left(kmer,
                                                     bundle,
                                                     new_kmers,
                                                     induced);
        }
        return 0;
    }

    uint8_t _find_induced_decision_nodes_left(kmer_t kmer,
                                              NeighborBundle& neighbors,
                                              set<hash_t>& new_kmers,
                                              deque<DecisionKmer>& induced) {

        // decision k-mers which are also new k-mers
        // cannot split existing unitigs. however, they can induce
        // an existing k-mer to be a decision k-mer, which can split
        // existing unitigs. so, we filter out neighbors of
        // the new decision k-mer which already exist and are already known
        // to the cDBG to be decision k-mers

        pdebug("Attempt left d-node induction from " << kmer);

        uint8_t n_found = 0;
        for (auto lneighbor : neighbors.first) {
            if (new_kmers.count(lneighbor.hash) ||
                cdbg->has_dnode(lneighbor.hash)) {
                continue;
            }
            NeighborBundle inductee_neighbors;
            if (get_decision_neighbors(lneighbor.kmer,
                                       inductee_neighbors,
                                       new_kmers)) {

                pdebug("Found induced d-node: " << lneighbor.hash << ", " << lneighbor.kmer);
                induced.push_back(make_pair(lneighbor, inductee_neighbors));
                ++n_found;
                //_build_dnode(lneighbor);
                //_split_unode(lneighbor, neighbors, neighbor_mask); 
            }
            inductee_neighbors.first.clear();
            inductee_neighbors.second.clear();
        }

        return n_found;
    }

    uint8_t _find_induced_decision_nodes_right(kmer_t kmer,
                                         set<hash_t>& new_kmers,
                                         deque<DecisionKmer>& induced) {

        pdebug("Prepare to attempt right induction on " << kmer);
        this->set_cursor(kmer.kmer);
        NeighborBundle bundle;
        bundle.second = find_right_kmers();
        
        if (bundle.second.size()) {
            return _find_induced_decision_nodes_right(kmer,
                                                      bundle,
                                                      new_kmers,
                                                      induced);
        }
        return 0;
    }

    uint8_t _find_induced_decision_nodes_right(kmer_t kmer,
                                         NeighborBundle& neighbors,
                                         set<hash_t>& new_kmers,
                                         deque<DecisionKmer>& induced) {

        // see _induce_decision_nodes_left for information

        pdebug("Attempt right d-node induction from " << kmer.kmer
                << ", " << kmer.hash);

        uint8_t n_found = 0;
        for (auto rneighbor : neighbors.second) {
            if (new_kmers.count(rneighbor.hash) ||
                cdbg->has_dnode(rneighbor.hash)) {
                continue;
            }

            NeighborBundle inductee_neighbors;
            if (get_decision_neighbors(rneighbor.kmer,
                                       inductee_neighbors,
                                       new_kmers)) {
                // induced decision k-mer
                pdebug("Found induced d-node: " << rneighbor.hash << ", " << rneighbor.kmer);
                induced.push_back(make_pair(rneighbor, inductee_neighbors));
                ++n_found;
                //_build_dnode(rneighbor);
                //_split_unode(rneighbor, neighbors, neighbor_mask); 
            }
            inductee_neighbors.first.clear();
            inductee_neighbors.second.clear();
        }

        return n_found;
    }

    virtual void _build_dnode(kmer_t kmer) {
        cdbg->build_dnode(kmer.hash, kmer.kmer);
    }

    uint8_t _add_neighbor_bundle(NeighborBundle& bundle) {
        return 0;    
    }
};

/*
template <class GraphType>
class AsyncStreamingCompactor : public StreamingCompactor<GraphType>,
                                public EventNotifier {

public:

    AsyncCDBG * acdbg;

    AsyncStreamingCompactor(GraphType * dbg,
                            uint64_t minimizer_window_size=8)
        : EventNotifier()
    {
        register_listener(static_cast<EventListener*>(acdbg));
        acdbg = new AsyncCDBG(dbg->K());
        StreamingCompactor<GraphType>::StreamingCompactor(dbg, minimizer_window_size, acdbg);
    }

    ~AsyncStreamingCompactor() {
        _cerr("AsyncStreamingCompactor: waiting for cDBG to finish updating.");
        //cdbg->wait_on_processing(0);

        delete acdbg;
    }

    void wait_on_updates() {
        acdbg->wait_on_processing(0);
    }

    void notify_build_dnode(hash_t hash, const string& kmer) {
        auto event = make_shared<BuildDNodeEvent>();
        event->hash = hash;
        event->kmer = kmer;
        this->notify(event);
    }

    void notify_build_unode(const string& sequence,
                            HashVector& tags,
                            hash_t left_end,
                            hash_t right_end) {
        auto event = make_shared<BuildUNodeEvent>();
        event->tags = tags;
        event->sequence = sequence;
        event->left_end = left_end;
        event->right_end = right_end;
        this->notify(event);
    }

    void notify_delete_unode(id_t node_id) {
        auto event = make_shared<DeleteUNodeEvent>();
        event->node_id = node_id;
        this->notify(event);
    }
};
*/
}
}
#undef pdebug
#endif
