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


#include "boink/assembly.hh"
#include "boink/hashing.hh"
#include "boink/dbg.hh"
#include "boink/cdbg.hh"
#include "boink/minimizers.hh"
#include "boink/event_types.hh"


namespace boink {
using namespace boink::event_types;


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

struct StreamingCompactorReport {
    uint64_t n_full;
    uint64_t n_tips;
    uint64_t n_islands;
    uint64_t n_unknown;
    uint64_t n_trivial;
    uint64_t n_dnodes;
    uint64_t n_unodes;
    uint64_t n_updates;
    uint64_t n_tags;
    uint64_t n_unique;
    double   estimated_fp;
};


struct compact_segment {
    hash_t left_anchor;
    hash_t right_anchor;
    bool has_left_decision;
    bool has_right_decision;
    string sequence;
    vector<hash_t> hashes;
};

template <class GraphType>
class StreamingCompactor : public AssemblerMixin<GraphType>,
                           public EventNotifier {

protected:

    uint64_t _minimizer_window_size;

public:

    using ShifterType = typename GraphType::shifter_type;
    using AssemblerType = AssemblerMixin<GraphType>;
    using AssemblerType::seen;
    using AssemblerType::get_left;
    using AssemblerType::get_right;
    using AssemblerType::degree_left;
    using AssemblerType::degree_right;
    using AssemblerType::count_nodes;
    using AssemblerType::filter_nodes;

    GraphType * dbg;
    cDBG cdbg;

    StreamingCompactor(GraphType * dbg,
                       uint64_t minimizer_window_size=8)
        : AssemblerMixin<GraphType>(dbg),
          EventNotifier(),
          _minimizer_window_size(minimizer_window_size),
          dbg(dbg),
          cdbg(dbg->K())
    {
        register_listener(static_cast<EventListener*>(&cdbg));
    }

    ~StreamingCompactor() {
        _cerr("StreamingCompactor: waiting for cDBG to finish updating.");
        cdbg.wait_on_processing(0);

        // make sure nothing else has a lock on the cdbg
        auto dlock = cdbg.lock_dnodes();
        auto ulock = cdbg.lock_unodes();
    }

    void wait_on_updates() {
        cdbg.wait_on_processing(0);
    }
    
    StreamingCompactorReport* get_report() {
        StreamingCompactorReport * report = new StreamingCompactorReport();
        report->n_full = cdbg.meta_counter.full_count;
        report->n_tips = cdbg.meta_counter.tip_count;
        report->n_islands = cdbg.meta_counter.island_count;
        report->n_unknown = cdbg.meta_counter.unknown_count;
        report->n_trivial = cdbg.meta_counter.trivial_count;
        report->n_dnodes = cdbg.n_decision_nodes();
        report->n_unodes = cdbg.n_unitig_nodes();
        report->n_tags = cdbg.n_tags();
        report->n_updates = cdbg.n_updates();
        report->n_unique = dbg->n_unique();
        report->estimated_fp = dbg->estimated_fp();

        return report;
    }

    string compactify(const string& seed) {
        Path path;
        this->set_cursor(seed);
        this->get_cursor(path);
        compactify_left(path);
        this->set_cursor(seed);
        compactify_right(path);

        return this->to_string(path);
    }

    void compactify_right(Path& path, InteriorMinimizer<hash_t>& minimizer) {
        this->seen.clear();
        this->seen.insert(this->get());
        minimizer.update(this->get());
        if (degree_left() > 1) return; // make sure we don't start at dnode
        
        shift_t next;
        while (get_right(next) // because get_right implicit checks rdegree
               && !this->seen.count(next.hash)) {
            path.push_back(next.symbol);
            this->seen.insert(next.hash);
            minimizer.update(next.hash);
            if (degree_left() > 1) break;
        }
    }

    void compactify_right(Path& path) {
        this->seen.clear();
        this->seen.insert(this->get());
        if (degree_left() > 1) return; // make sure we don't start at dnode
        
        shift_t next;
        while (get_right(next) // because get_right implicit checks rdegree
               && !this->seen.count(next.hash)) {
            path.push_back(next.symbol);
            this->seen.insert(next.hash);
            if (degree_left() > 1) break;
        }
    }

    void compactify_left(Path& path, InteriorMinimizer<hash_t> minimizer) {
        this->seen.clear();
        this->seen.insert(this->get());
        minimizer.update(this->get());
        if (degree_right() > 1) return; // don't start at dnode

        shift_t next;
        while (get_left(next) // because get_left implicitly checks ldegree
               && !this->seen.count(next.hash)) {
            path.push_front(next.symbol);
            this->seen.insert(next.hash);
            minimizer.update(next.hash);
            if (degree_right() > 1) break;
        }
    }

    void compactify_left(Path& path) {
        this->seen.clear();
        this->seen.insert(this->get());
        if (degree_right() > 1) return; // don't start at dnode

        shift_t next;
        while (get_left(next) // because get_left implicitly checks ldegree
               && !this->seen.count(next.hash)) {
            path.push_front(next.symbol);
            this->seen.insert(next.hash);
            if (degree_right() > 1) break;
        }
    }

    void notify_build_dnode(hash_t hash, const string& kmer) {
        auto event = make_shared<BuildDNodeEvent>();
        event->hash = hash;
        event->kmer = kmer;
        this->notify(event);
    }

    void notify_build_unode(const string& sequence,
                            HashVector& tags,
                            bool has_left,
                            hash_t left_end,
                            bool has_right,
                            hash_t right_end) {
        auto event = make_shared<BuildUNodeEvent>();
        event->tags = tags;
        event->sequence = sequence;
        event->has_left = has_left;
        event->left_end = left_end;
        event->has_right = has_right;
        event->right_end = right_end;
        this->notify(event);
    }

    void notify_delete_unode(id_t node_id) {
        auto event = make_shared<DeleteUNodeEvent>();
        event->node_id = node_id;
        this->notify(event);
    }

    bool is_decision_kmer(const string& node,
                          uint8_t& degree) {
        this->set_cursor(node);
        return is_decision_kmer(degree);
    }

    bool is_decision_kmer(const string& node) {
        this->set_cursor(node);
        return this->degree_left() > 1 || this->degree_right() > 1;
    }

    bool is_decision_kmer(uint8_t& degree) {
        uint8_t ldegree, rdegree;
        ldegree = this->degree_left();
        rdegree = this->degree_right();
        degree = ldegree + rdegree;
        return ldegree > 1 || rdegree > 1;
    }

    void find_decision_kmers(const string& sequence,
                             vector<uint32_t>& decision_positions,
                             HashVector& decision_hashes,
                             vector<NeighborBundle>& decision_neighbors) {

        KmerIterator<typename GraphType::shifter_type> iter(sequence, this->_K);
        size_t pos = 0;
        while(!iter.done()) {
            hash_t h = iter.next();
            NeighborBundle neighbors;
            if (get_decision_neighbors(iter.shifter,
                                       sequence.substr(pos, this->_K),
                                       neighbors)) {

                decision_neighbors.push_back(neighbors);
                decision_positions.push_back(pos);
                decision_hashes.push_back(h);
            }

            ++pos;
        }
    }

    void update_sequence(const string& sequence) {

    }

    bool get_decision_neighbors(typename GraphType::shifter_type& shifter,
                                const string& root_kmer,
                                NeighborBundle& result) {

        vector<shift_t> left_neighbors = filter_nodes(shifter.gather_left());
        vector<shift_t> right_neighbors = filter_nodes(shifter.gather_right());

        if (left_neighbors.size() > 1 || right_neighbors.size() > 1) {
            KmerVector left_kmers, right_kmers;
            for (auto neighbor : left_neighbors) {
                left_kmers.push_back(kmer_t(neighbor.hash,
                                            neighbor.symbol
                                            + root_kmer.substr(0, this->_K-1)));
            }
            for (auto neighbor : right_neighbors) {
                right_kmers.push_back(kmer_t(neighbor.hash,
                                             root_kmer.substr(1, this->_K-1)
                                             + neighbor.symbol));
                                      
            }
            result = std::make_pair(left_kmers, right_kmers);
            return true;
        } else {
            return false;
        }
    }

    //void find_new_segments(const string& sequence,
    //                       vector<vector<hash_t>>& new_segment_hashes,
    //                       vector<string>& new_segment_sequences,
    //                       deque<kmer_t>& decision_kmers,
    //                       deque<NeighborBundle>& decision_neighbors
    //                       ) {


    compact_segment init_segment(hash_t left_anchor,
                                 bool has_left_decision) {
        compact_segment segment;
        segment.left_anchor = left_anchor;
        segment.has_left_decision = has_left_decision;
        segment.has_right_decision = false;

        return segment;
    }

    void finish_segment(const string& sequence,
                        compact_segment& segment,
                        size_t start,
                        size_t end,
                        hash_t right_anchor,
                        deque<compact_segment>& segments) {
        segment.sequence = sequence.substr(start, end - start + this->_K);
        segment.right_anchor = right_anchor;
        segments.push_back(segment);
        pdebug("segment ended, start=" << start
               << " cur_pos=" << end << ", seq=" << segment.sequence);
    }
    
    void find_new_segments(const string& sequence,
                           deque<compact_segment>& segments,
                           deque<NeighborBundle>& decision_neighbors
                           ) {


        vector<hash_t> kmer_hashes;
        vector<bool> kmer_new;
        dbg->add_sequence(sequence, kmer_hashes, kmer_new);
        pdebug("k-mers new: " << kmer_new);
        
        size_t pos = 0;
        size_t segment_start_pos = 0;
        hash_t prev_hash = kmer_hashes.front(), cur_hash;
        bool cur_new = false, prev_new = false, is_decision = false;
        string kmer_seq;

        compact_segment current_segment;
        typename GraphType::shifter_type shifter(sequence, this->_K);
        for (cur_hash : kmer_hashes) {
            cur_new = kmer_new[pos];

            if(cur_new) {
                kmer_seq = sequence.substr(pos, this->_K);

                if(!prev_new) {
                    pdebug("old -> new (pos=" << pos << ")");
                    current_segment = init_segment(prev_hash, false);
                    segment_start_pos = pos;
                    shifter.set_cursor(kmer_seq);
                } else {
                    shifter.shift_right(kmer_seq.back());
                }
                
                NeighborBundle neighbors;
                neighbors.first.clear();
                neighbors.second.clear();
                if (get_decision_neighbors(shifter,
                                           kmer_seq,
                                           neighbors)) {
                    decision_neighbors.push_back(neighbors);
                }
                is_decision = neighbors.first.size() > 1 || neighbors.second.size() > 1;
                if(is_decision) {
                    pdebug("new k-mer & decision" << shifter.get_cursor()
                           << ", " << kmer_seq);
                    current_segment.has_right_decision = true;
                    if (pos > 0) {
                        // don't want to produce a null segment if the decision node
                        // is the first in the sequence
                        finish_segment(sequence, current_segment,
                                       segment_start_pos, pos, cur_hash,
                                       segments);
                    }
                    current_segment = init_segment(cur_hash, true);
                    segment_start_pos = pos;
                }
                current_segment.hashes.push_back(cur_hash);
            } else if (prev_new) {
                pdebug("new -> old");
                finish_segment(sequence, current_segment, 
                               segment_start_pos, pos - 1, cur_hash,
                               segments);
            } 

            ++pos;
            prev_hash = cur_hash;
            prev_new = cur_new;
        }

        if (cur_new) {
            pdebug("sequence ended");
            current_segment.has_right_decision = is_decision;
            finish_segment(sequence, current_segment,
                           segment_start_pos, pos, cur_hash,
                           segments);
        }

        //pdebug("Segments: " << new_segment_sequences);
    }

    void update_from_segments(const string& sequence,
                              vector<vector<hash_t>>& new_segment_kmers,
                              vector<string>& new_segment_sequences,
                              deque<kmer_t>&& decision_kmers,
                              deque<NeighborBundle>& decision_neighbors
                              ) {

        size_t dnode_counter = 0;
        for (auto segment_kmers : new_segment_kmers) {
            // handle front k-mer
            if (segment_kmers.front() !=
                decision_kmers.front().hash) {

                connect_segment_end(segment_kmers.front(), DIR_LEFT);
            } else {
                new_decision_kmer(decision_kmers.front(),
                                  decision_neighbors.front());
                decision_kmers.pop_front();
                decision_neighbors.pop_front();
            }

            // handle middle k-mers
            // only k-mers we do anything with are decision k-mers
            if (segment_kmers.size() > 2) {
                auto kmer_iter = segment_kmers.begin() + 1;
                while(kmer_iter != segment_kmers.end() - 1) {
                    auto kmer = *kmer_iter;

                    if (kmer == decision_kmers.front().hash) {
                        new_decision_kmer(decision_kmers.front(),
                                          decision_neighbors.front());
                        decision_kmers.pop_front();
                        decision_neighbors.pop_front();
                    }
                }
            }

            // handle back k-mer
        }
    }

    void new_decision_kmer() {
        /* Handle a new k-mer that is also
         * a decision node.
         */
    }

    void induce_decision_node() {
        /* Handle an induced decision node (an existing
         * k-mer converted to a d-node), which splits an
         * existing unitig.
         */
    }

    void connect_segment_end(kmer_t kmer,
                              direction_t direction) {
        
    }
};

}
#undef pdebug
#endif
