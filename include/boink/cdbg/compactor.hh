/**
 * (c) Camille Scott, 2019
 * File   : compactor.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 18.07.2019
 */

#ifndef BOINK_COMPACTOR_HH
#define BOINK_COMPACTOR_HH

#include <assert.h>
#include <cstdint>

#include "boink/traversal.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/dbg.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/nibblestorage.hh"
#include "boink/cdbg/cdbg.hh"
#include "boink/minimizers.hh"

#include "boink/event_types.hh"

#include "boink/reporting/reporters.hh"
#include "boink/processors.hh"


namespace boink {
namespace cdbg {


# ifdef DEBUG_CPTR
#   define pdebug(x) do { std::cerr << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif


template <class GraphType>
struct StreamingCompactor {

    using ShifterType    = typename GraphType::shifter_type;
    using TraversalType  = Traverse<GraphType>;
    using cDBGType       = typename cDBG<GraphType>::Graph;
    using CompactNode    = typename cDBG<GraphType>::CompactNode;
    using UnitigNode     = typename cDBG<GraphType>::UnitigNode;
    using DecisionNode   = typename cDBG<GraphType>::DecisionNode;
    using State          = typename TraversalType::State;
    using NeighborBundle = typename TraversalType::NeighborBundle;

    typedef GraphType                         graph_type;

    typedef ShifterType                       shifter_type;
    typedef typename shifter_type::hash_type  hash_type;
    typedef typename shifter_type::kmer_type  kmer_type;
    typedef typename shifter_type::shift_type shift_type;

    typedef std::pair<kmer_type, NeighborBundle> DecisionKmer;


    /* Represents a segment of new k-mers from a sequence, relative to the current
     * state of the cDBG. Can either be: a null segment representing a portion
     * of the sequence already in the graph; a decision k-mer; or a new unitig or
     * portion of a unitig.
     */

    struct compact_segment {
        hash_type left_anchor;
        hash_type right_anchor;
        hash_type left_flank;
        hash_type right_flank;

        // whether this segment represents a decision k-mer
        bool is_decision_kmer;
        // start position of the segment within the originating sequence
        size_t start_pos;
        // length of the segment sequence (from beginning of first k-mer to
        // end of last k-mer)
        size_t length;
        // tags associated with this segment
        std::vector<hash_type> tags;

        // the default constructor creates a null segment
        compact_segment()
            : left_anchor(0),
              right_anchor(0),
              left_flank(0),
              right_flank(0),
              is_decision_kmer(false),
              start_pos(0),
              length(0) 
        {
        }

        compact_segment(hash_type left_anchor,
                        hash_type right_anchor,
                        bool is_decision_kmer,
                        size_t start_pos,
                        size_t length)
            : left_anchor(left_anchor),
              right_anchor(right_anchor),
              is_decision_kmer(is_decision_kmer),
              start_pos(start_pos),
              length(length) 
        {
        }

        compact_segment(const compact_segment &obj)
            : left_anchor(obj.left_anchor),
              right_anchor(obj.right_anchor),
              left_flank(obj.left_flank),
              right_flank(obj.right_flank),
              is_decision_kmer(obj.is_decision_kmer),
              start_pos(obj.start_pos),
              length(obj.length)
        {
        }

        // return if the segment is default constructed / null
        // used as a delimiter token between connected runs of segments
        const bool is_null() const {
            return length == 0 && (left_anchor == right_anchor) && !is_decision_kmer;
        }

        friend inline std::ostream& operator<<(std::ostream& os, const compact_segment& segment) {
            os << "<compact_segment left_flank=" << segment.left_flank
               << " left_anchor=" << segment.left_anchor
               << " right_anchor=" << segment.right_anchor
               << " right_flank=" << segment.right_flank
               << " start=" << segment.start_pos
               << " length=" << segment.length
               << ">";
            return os;
        }
    };


    struct Report {
        uint64_t n_full;
        uint64_t n_tips;
        uint64_t n_islands;
        uint64_t n_trivial;
        uint64_t n_circular;
        uint64_t n_loops;
        uint64_t n_dnodes;
        uint64_t n_unodes;
        uint64_t n_updates;
        uint64_t n_splits;
        uint64_t n_merges;
        uint64_t n_extends;
        uint64_t n_clips;
        uint64_t n_deletes;
        uint64_t n_circular_merges;
        uint64_t n_tags;
        uint64_t n_unique;
        double   estimated_fp;
    };

    class Compactor : public TraversalType::dBG,
                      public events::EventNotifier {

    protected:

        uint64_t _minimizer_window_size;

    public:



        using TraversalType::dBG::filter_nodes;
        using TraversalType::dBG::find_left_kmers;
        using TraversalType::dBG::find_right_kmers;
        using TraversalType::dBG::gather_left;
        using TraversalType::dBG::gather_right;
        using TraversalType::dBG::traverse_left;
        using TraversalType::dBG::traverse_right;
        using TraversalType::dBG::get_decision_neighbors;

        typedef ShifterType   shifter_type;
        typedef GraphType     graph_type;
        typedef cDBGType      cdbg_type;

        shared_ptr<GraphType> dbg;
        shared_ptr<cDBGType> cdbg;

        Compactor(shared_ptr<GraphType> dbg,
                  uint64_t minimizer_window_size=8)
            : TraversalType::dBG(dbg->K()),
              EventNotifier(),
              _minimizer_window_size(minimizer_window_size),
              dbg(dbg)
        {
            this->cdbg = make_shared<cDBGType>(dbg,
                                               minimizer_window_size);
        }

        ~Compactor() {
            //cdbg->wait_on_processing(0);

            // make sure nothing else has a lock on the cdbg
            auto lock = cdbg->lock_nodes();
        }

        static std::shared_ptr<Compactor> build(std::shared_ptr<GraphType> dbg,
                                                uint64_t minimizer_window_size=8) {
            return std::make_shared<Compactor>(dbg, minimizer_window_size);
        }
        
        Report get_report() {
            Report report;
            report.n_full            = cdbg->metrics->n_full;
            report.n_tips            = cdbg->metrics->n_tips;
            report.n_islands         = cdbg->metrics->n_islands;
            report.n_trivial         = cdbg->metrics->n_trivial;
            report.n_circular        = cdbg->metrics->n_circular;
            report.n_loops           = cdbg->metrics->n_loops;
            report.n_dnodes          = cdbg->metrics->n_dnodes;
            report.n_unodes          = cdbg->metrics->n_unodes;
            report.n_tags            = cdbg->n_tags();
            report.n_updates         = cdbg->n_updates();
            report.n_splits          = cdbg->metrics->n_splits;
            report.n_merges          = cdbg->metrics->n_merges;
            report.n_extends         = cdbg->metrics->n_extends;
            report.n_clips           = cdbg->metrics->n_clips;
            report.n_deletes         = cdbg->metrics->n_deletes;
            report.n_circular_merges = cdbg->metrics->n_circular_merges;
            report.n_unique          = dbg->n_unique();
            //report.estimated_fp      = dbg->estimated_fp();
            report.estimated_fp      = 0;

            return report;
        }

        size_t insert_sequence(const std::string& sequence,
                               std::shared_ptr<std::vector<hash_type>> hashes = nullptr) {
            std::set<hash_type> new_kmers;
            std::deque<compact_segment> segments;
            std::set<hash_type> new_decision_kmers;
            std::deque<NeighborBundle> decision_neighbors;
            //std::vector<hash_type> hashes;
            if (hashes == nullptr) {
                hashes = std::make_shared<std::vector<hash_type>>();
            }

            find_new_segments(sequence,
                              *hashes,
                              new_kmers,
                              segments,
                              new_decision_kmers,
                              decision_neighbors);

            update_from_segments(sequence,
                                 new_kmers,
                                 segments,
                                 new_decision_kmers,
                                 decision_neighbors);

            for (auto& h : *hashes) {
                dbg->insert(h);
            }

            return hashes->size();
        }

        compact_segment init_segment(hash_type left_anchor,
                                     hash_type left_flank,
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
                            hash_type right_anchor,
                            hash_type right_flank,
                            std::deque<compact_segment>& segments) {
            
            segment.length = end - segment.start_pos + this->_K;
            segment.right_anchor = right_anchor;
            segment.right_flank = right_flank;
            segments.push_back(segment);
            pdebug("Finished segment: " << segment);
        }

        void finish_decision_segment(compact_segment& segment,
                                     std::deque<compact_segment>& segments) {
            segment.length = this->_K;
            segment.right_anchor = segment.left_anchor;
            segment.right_flank = segment.left_flank;
            segment.is_decision_kmer = true;
            segments.push_back(segment);

            pdebug("Finished decision segment: " << segment);       
        }

        std::deque<compact_segment> find_new_segments(const std::string& sequence) {
        
            std::deque<compact_segment> result;
            find_new_segments(sequence, result);
            return result;
        }

        void find_new_segments(const std::string& sequence,
                               std::deque<compact_segment>& result) {
            // convenience function for cython land
            std::vector<hash_type> hashes;
            std::set<hash_type> new_kmers;
            std::set<hash_type> new_decision_kmers;
            std::deque<NeighborBundle> decision_neighbors;

            find_new_segments(sequence,
                              hashes,
                              new_kmers,
                              result,
                              new_decision_kmers,
                              decision_neighbors);

        }

        void find_new_segments(const std::string& sequence,
                               std::vector<hash_type>& hashes,
                               std::set<hash_type>& new_kmers,
                               std::deque<compact_segment>& segments,
                               std::set<hash_type>& new_decision_kmers,
                               std::deque<NeighborBundle>& decision_neighbors) {

            pdebug("FIND SEGMENTS: " << sequence);

            hashing::KmerIterator<ShifterType> kmers(sequence, this->_K);
            hash_type prev_hash, cur_hash;
            size_t pos = 0;
            bool cur_new = false, prev_new = false, cur_seen = false, prev_seen = false;

            std::deque<compact_segment> preprocess;
#ifdef DEBUG_CPTR
            std::vector<count_t> counts;
#endif
            compact_segment current_segment; // start null
            std::vector<typename TraversalType::dBG> segment_shifters;
            while(!kmers.done()) {
                cur_hash = kmers.next();
                cur_new = this->dbg->query(cur_hash) == 0;
                cur_seen = new_kmers.count(cur_hash);
                hashes.push_back(cur_hash);
#ifdef DEBUG_CPTR
                counts.push_back(this->dbg->query(cur_hash));
#endif

                if(cur_new && !cur_seen) {

                    if(!prev_new || prev_seen) {
                        pdebug("old -> new (pos=" << pos << ")");

                        preprocess.push_back(current_segment);
                        current_segment = init_segment(cur_hash,
                                                       prev_hash,
                                                       pos);
                        segment_shifters.emplace_back(typename TraversalType::dBG(*(kmers.shifter)));
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
                hash_type right_flank = cur_hash;
                std::vector<shift_type> rneighbors = filter_nodes(dbg.get(),
                                                               kmers.shifter->gather_right(),
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
                std::vector<shift_type> lneighbors = filter_nodes(dbg.get(),
                                                               segment_shifters[0].gather_left(),
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

                std::deque<compact_segment> decision_segments;
                size_t pos = segment.start_pos;
                size_t suffix_pos = pos + this->_K - 1;
                while (1) {
                    NeighborBundle neighbors;
                    if (shifter_iter->get_decision_neighbors(dbg.get(),
                                                             neighbors,
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

       void update_from_segments(const std::string& sequence,
                                  std::set<hash_type>& new_kmers,
                                  std::deque<compact_segment>& segments,
                                  std::set<hash_type> & new_decision_kmers,
                                  std::deque<NeighborBundle>& decision_neighbors
                                  ) {

            if (segments.size() < 3) {
                pdebug("No segments.");
                return;
            }

            pdebug(segments.size() << " segments.");

            // First find all induced decision k-kmers
            // Have to wait until all are found to proceed with induction process
            // in the cDBG though
            std::deque<DecisionKmer> induced;
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
                    kmer_type decision_kmer(v.left_anchor,
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
                        kmer_type root(v.left_anchor,
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
                        kmer_type root(v.right_anchor,
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

        void _induce_decision_nodes(std::deque<DecisionKmer>& induced_decision_kmers,
                                    std::set<hash_type>& new_kmers) {

            std::set<hash_type> induced_decision_kmer_hashes;

            pdebug("Perform induction on " << induced_decision_kmers.size() <<
                   " new decision k-mers");
            for (auto dkmer : induced_decision_kmers) {
                _build_dnode(dkmer.first);
                induced_decision_kmer_hashes.insert(dkmer.first.hash);
            }

            std::set<hash_type> processed;
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

        virtual bool _try_split_unode(kmer_type root,
                                  NeighborBundle& neighbors,
                                  std::set<hash_type>& new_kmers,
                                  std::set<hash_type>& induced_decision_kmer_hashes,
                                  std::set<hash_type>& processed) {
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

                hash_type new_end;
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

            std::vector<kmer_type> lfiltered;
            std::copy_if(neighbors.first.begin(),
                         neighbors.first.end(),
                         std::back_inserter(lfiltered),
                         [&] (kmer_type neighbor) { return
                            !new_kmers.count(neighbor.hash) &&
                            !processed.count(neighbor.hash);
                         });

            std::vector<kmer_type> rfiltered;
            std::copy_if(neighbors.second.begin(),
                         neighbors.second.end(),
                         std::back_inserter(rfiltered),
                         [&] (kmer_type neighbor) { return
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
                Path path;
                auto stop_state = this->traverse_left(dbg.get(), start.kmer, path, processed);
                hash_type end_hash = stop_state.second;

                if (stop_state.first != State::STOP_SEEN) {

                    if (stop_state.first == State::DECISION_FWD) {

                        shift_type right;
                        auto n_found = this->try_traverse_right(dbg.get(), right);
                        if (n_found == 1) {
                            end_hash = right.hash;
                            path.pop_front();
                        }
                    }
                    unode_to_split = cdbg->query_unode_end(end_hash);
                    size_t split_point = path.size() - this->_K  + 1;
                    hash_type left_unode_new_right = start.hash;

                    pdebug("split point is " << split_point <<
                            " new_right is " << left_unode_new_right
                           << " root was " << root.hash);

                    hash_type right_unode_new_left = this->hash(unode_to_split->sequence.c_str() + 
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
                std::cout << "rstart: " << start.kmer << std::endl;
                Path path;
                auto stop_state = this->traverse_right(dbg.get(), start.kmer, path, processed);
                hash_type end_hash = stop_state.second;

                if (stop_state.first != State::STOP_SEEN) {

                    if (stop_state.first == State::DECISION_FWD) {

                        shift_type right;
                        auto n_found = this->try_traverse_left(dbg.get(), right);
                        if (n_found == 1) {
                            end_hash = right.hash;
                            path.pop_back();
                        }
                    }

                    unode_to_split = cdbg->query_unode_end(end_hash);
                    size_t split_point = unode_to_split->sequence.size()
                                                         - path.size()
                                                         - 1;
                    hash_type new_right = this->hash(unode_to_split->sequence.c_str() + 
                                                  split_point - 1);
                    hash_type new_left = start.hash;
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
                                   const std::string& sequence) {

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
                std::string trimmed_seq;
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

        uint8_t _find_induced_decision_nodes(kmer_type kmer,
                                          NeighborBundle& neighbors,
                                          std::set<hash_type>& new_kmers,
                                          std::deque<DecisionKmer>& induced) {

            return _find_induced_decision_nodes_left(kmer, neighbors, new_kmers, induced) +
                   _find_induced_decision_nodes_right(kmer, neighbors, new_kmers, induced);
        }

        uint8_t _find_induced_decision_nodes_left(kmer_type kmer,
                                                  std::set<hash_type>& new_kmers,
                                                  std::deque<DecisionKmer>& induced) {

            pdebug("Prepare to attempt left induction on " << kmer);
            this->set_cursor(kmer.kmer);
            NeighborBundle bundle;
            bundle.first = find_left_kmers(dbg.get());

            if (bundle.first.size()) {
                return _find_induced_decision_nodes_left(kmer,
                                                         bundle,
                                                         new_kmers,
                                                         induced);
            }
            return 0;
        }

        uint8_t _find_induced_decision_nodes_left(kmer_type kmer,
                                                  NeighborBundle& neighbors,
                                                  std::set<hash_type>& new_kmers,
                                                  std::deque<DecisionKmer>& induced) {

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
                if (get_decision_neighbors(dbg.get(),
                                           lneighbor.kmer,
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

        uint8_t _find_induced_decision_nodes_right(kmer_type kmer,
                                             std::set<hash_type>& new_kmers,
                                             std::deque<DecisionKmer>& induced) {

            pdebug("Prepare to attempt right induction on " << kmer);
            this->set_cursor(kmer.kmer);
            NeighborBundle bundle;
            bundle.second = find_right_kmers(dbg.get());
            
            if (bundle.second.size()) {
                return _find_induced_decision_nodes_right(kmer,
                                                          bundle,
                                                          new_kmers,
                                                          induced);
            }
            return 0;
        }

        uint8_t _find_induced_decision_nodes_right(kmer_type kmer,
                                             NeighborBundle& neighbors,
                                             std::set<hash_type>& new_kmers,
                                             std::deque<DecisionKmer>& induced) {

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
                if (get_decision_neighbors(dbg.get(),
                                           rneighbor.kmer,
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

        virtual void _build_dnode(kmer_type kmer) {
            cdbg->build_dnode(kmer.hash, kmer.kmer);
        }

        uint8_t _add_neighbor_bundle(NeighborBundle& bundle) {
            return 0;    
        }

        void reverse_complement_cdbg() {
            for (auto it = cdbg->dnodes_begin(); it != cdbg->dnodes_end(); ++it) {
                auto rc_sequence = it->second->revcomp();
                insert_sequence(rc_sequence);
            }
            for (auto it = cdbg->unodes_begin(); it != cdbg->unodes_end(); ++it) {
                auto rc_sequence = it->second->revcomp();
                insert_sequence(rc_sequence);
            }

            std::vector<CompactNode*>  rc_nodes_to_delete;
            spp::sparse_hash_set<id_t> visited;
            for (auto it = cdbg->unodes_begin(); it != cdbg->unodes_end(); ++it) {
                if (visited.count(it->second->node_id)) {
                    continue;
                }
                auto rc_node = cdbg->find_rc_cnode(it->second.get());
                rc_nodes_to_delete.push_back(rc_node);
                
                visited.insert(rc_node->node_id);
                visited.insert(it->second->node_id);
            }

            for (auto it = cdbg->dnodes_begin(); it != cdbg->dnodes_end(); ++it) {
                if (visited.count(it->second->node_id)) {
                    continue;
                }
                auto rc_node = cdbg->find_rc_cnode(it->second.get());
                if (rc_node->node_id != it->second->node_id) {
                    rc_nodes_to_delete.push_back(rc_node);
                }

                visited.insert(rc_node->node_id);
                visited.insert(it->second->node_id);
            }
         
            for (auto node : rc_nodes_to_delete) {
                if (node->meta() == DECISION) {
                    cdbg->delete_dnode((DecisionNode*)node);
                } else {
                    cdbg->delete_unode((UnitigNode*)node);
                }
            }
        }
    };


    using Processor = InserterProcessor<Compactor>;


    class Reporter: public reporting::SingleFileReporter {

    protected:

        std::shared_ptr<Compactor> compactor;

    public:

        Reporter(std::shared_ptr<Compactor> compactor,
                 const std::string& output_filename)
            : SingleFileReporter(output_filename, "StreamingCompactor::Reporter"),
              compactor(compactor)
        {    
            _cerr(this->THREAD_NAME << " reporting at FINE interval.");

            this->msg_type_whitelist.insert(events::MSG_TIME_INTERVAL);

            _output_stream << "read_n,n_full,n_tips,n_islands,n_trivial"
                              ",n_circular,n_loops,n_dnodes,n_unodes,n_tags,"
                              "n_updates,n_splits,n_merges,n_extends,n_clips,"
                              "n_deletes,n_circular_merges,n_unique,estimated_fp" << std::endl;
        }

        static std::shared_ptr<Reporter> build(std::shared_ptr<Compactor> compactor,
                                               const std::string& output_filename) {
            return std::make_shared<Reporter>(compactor, output_filename);
        }

        virtual void handle_msg(shared_ptr<events::Event> event) {
            if (event->msg_type == events::MSG_TIME_INTERVAL) {
                auto _event = static_cast<events::TimeIntervalEvent*>(event.get());
                if (_event->level == events::TimeIntervalEvent::FINE ||
                    _event->level == events::TimeIntervalEvent::END) {
                    auto report = compactor->get_report();
                    _output_stream << _event->t << ","
                                   << report.n_full << ","
                                   << report.n_tips << ","
                                   << report.n_islands << ","
                                   << report.n_trivial << ","
                                   << report.n_circular << ","
                                   << report.n_loops << ","
                                   << report.n_dnodes << ","
                                   << report.n_unodes << ","
                                   << report.n_tags << ","
                                   << report.n_updates << ","
                                   << report.n_splits << ","
                                   << report.n_merges << ","
                                   << report.n_extends << ","
                                   << report.n_clips << ","
                                   << report.n_deletes << ","
                                   << report.n_circular_merges << ","
                                   << report.n_unique << ","
                                   << report.estimated_fp 
                                   << std::endl;
                }
            }
        }
    };

    template <class ParserType = parsing::FastxReader>
    class NormalizingCompactor : public FileProcessor<NormalizingCompactor<ParserType>,
                                                      ParserType> { 
    protected:


        bool median_count_at_least(const std::string&          sequence,
                                   unsigned int                cutoff,
                                   dBG<storage::ByteStorage,
                                       ShifterType>            * counts) {

            auto kmers = counts->get_hash_iter(sequence);
            unsigned int min_req = 0.5 + float(sequence.size() - counts->K() + 1) / 2;
            unsigned int num_cutoff_kmers = 0;

            // first loop:
            // accumulate at least min_req worth of counts before checking to see
            // if we have enough high-abundance k-mers to indicate success.
            for (unsigned int i = 0; i < min_req; ++i) {
                hash_type kmer = kmers->next();
                if (counts->query(kmer) >= cutoff) {
                    ++num_cutoff_kmers;
                }
            }


            // second loop: now check to see if we pass the threshold for each k-mer.
            if (num_cutoff_kmers >= min_req) {
                return true;
            }
            while(!kmers->done()) {
                hash_type kmer = kmers->next();
                if (counts->query(kmer) >= cutoff) {
                    ++num_cutoff_kmers;
                    if (num_cutoff_kmers >= min_req) {
                        return true;
                    }
                }
            }
            return false;
        }


        std::shared_ptr<Compactor>                              compactor;
        std::shared_ptr<GraphType>                              graph;

        std::unique_ptr<dBG<storage::ByteStorage, ShifterType>> counts;
        unsigned int                                            cutoff;
        size_t                                                  n_seq_updates;

        typedef FileProcessor<NormalizingCompactor<ParserType>,
                              ParserType> Base;

    public:

        using Base::process_sequence;
        using events::EventNotifier::register_listener;
        
        NormalizingCompactor(shared_ptr<Compactor> compactor,
                             unsigned int cutoff,
                             uint64_t fine_interval   = DEFAULT_INTERVALS::FINE,
                             uint64_t medium_interval = DEFAULT_INTERVALS::MEDIUM,
                             uint64_t coarse_interval = DEFAULT_INTERVALS::COARSE)
            : Base(fine_interval, medium_interval, coarse_interval),
              compactor(compactor),
              graph(compactor->dbg),
              cutoff(cutoff),
              n_seq_updates(0)
        {
            auto storage = storage::ByteStorage::build(100000000, 4);
            auto hasher = graph->get_hasher();
            counts = std::make_unique<dBG<storage::ByteStorage,
                                          ShifterType>>(hasher, storage);
        }

        static std::shared_ptr<NormalizingCompactor> build(std::shared_ptr<Compactor> compactor,
                                                           unsigned int cutoff,
                                                           uint64_t fine_interval   = DEFAULT_INTERVALS::FINE,
                                                           uint64_t medium_interval = DEFAULT_INTERVALS::MEDIUM,
                                                           uint64_t coarse_interval = DEFAULT_INTERVALS::COARSE) {
            return std::make_shared<NormalizingCompactor>(compactor, cutoff, fine_interval, medium_interval, coarse_interval);

        }

        void process_sequence(const parsing::Read& read) {

            if (median_count_at_least(read.cleaned_seq, cutoff, counts.get())) {
                return;
            }

            counts->insert_sequence(read.cleaned_seq);

            try {
                compactor->insert_sequence(read.cleaned_seq);
            } catch (hashing::InvalidCharacterException &e) {
                std::cerr << "WARNING: Bad sequence encountered at "
                          << this->_n_reads << ": "
                          << read.cleaned_seq << ", exception was "
                          << e.what() << std::endl;
                return;
            } catch (hashing::SequenceLengthException &e) {
                std::cerr << "NOTE: Skipped sequence that was too short: read "
                          << this->_n_reads << " with sequence "
                          << read.cleaned_seq 
                          << std::endl;
                return;
            } catch (std::exception &e) {
                std::cerr << "ERROR: Exception thrown at " << this->_n_reads 
                          << " with msg: " << e.what()
                          <<  std::endl;
                throw e;
            }

            ++n_seq_updates;
        }

        void report() {
            std::cerr << "\t" << n_seq_updates << " used for cDBG updates." << std::endl;
        }

    };

    
    class SolidCompactor : public events::EventNotifier {

    private:

        std::unique_ptr<dBG<storage::NibbleStorage,
                            ShifterType>> abund_filter;

    public:

        std::shared_ptr<Compactor>    compactor;
        std::shared_ptr<GraphType>    dbg;

        unsigned int                  min_abund;

        SolidCompactor(std::shared_ptr<Compactor> compactor,
                       unsigned int               min_abund,
                       uint64_t                   abund_table_size,
                       uint16_t                   n_abund_tables)
            : EventNotifier (),
              compactor     (compactor),
              dbg           (compactor->dbg),
              min_abund     (min_abund)
        {
            std::shared_ptr<storage::NibbleStorage> abund_storage = storage::NibbleStorage::build(abund_table_size,
                                                                                                  n_abund_tables);
            auto abund_hasher = dbg->get_hasher();
            abund_filter = std::make_unique<dBG<storage::NibbleStorage,
                                                ShifterType>>(abund_hasher,
                                                              abund_storage);
        }

        static std::shared_ptr<SolidCompactor> build(std::shared_ptr<Compactor> compactor,
                                                     unsigned int               min_abund,
                                                     uint64_t                   abund_table_size,
                                                     uint16_t                   n_abund_tables) {
            return std::make_shared<SolidCompactor>(compactor, min_abund, abund_table_size, n_abund_tables);
        }
                                                     

        std::vector<std::pair<size_t, size_t>> find_solid_segments(const std::string& sequence) {
            std::vector<hash_type>  hashes;
            std::vector<storage::count_t> counts;
            abund_filter->insert_sequence(sequence, hashes, counts);

            std::vector<std::pair<size_t, size_t>> segments;
            auto hashes_iter = hashes.begin();
            auto counts_iter = counts.begin();
            size_t start, end, pos = 0;
            hash_type cur_hash;
            bool   prev_solid = false, cur_solid = false;
            while(hashes_iter != hashes.end()) {
                cur_hash  = *hashes_iter;
                cur_solid = (*counts_iter) >= min_abund;
                if (cur_solid && !prev_solid) {
                    start = pos;
                } else if (prev_solid && !cur_solid) {
                    end = pos + compactor->K() - 1;
                    segments.emplace_back(start, end);
                }
                
                ++pos;
                ++hashes_iter;
                ++counts_iter;
                prev_solid = cur_solid;
            }

            if (cur_solid) {
                end = pos + compactor->K() - 1;
                segments.emplace_back(start, end);
            }

            return segments;
        }

        void insert_sequence(const std::string& sequence) {

            auto solid_segments = find_solid_segments(sequence);
            for (auto segment : solid_segments) {
                auto length = segment.second - segment.first;
                compactor->insert_sequence(sequence.substr(segment.first,
                                                           length));
            }
        }
              
                                    
    };

};


}
}
#undef pdebug
#endif
