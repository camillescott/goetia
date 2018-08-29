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

/* Packs up information on the current compaction state.
 */
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


/* Represents a segment of new k-mers from a sequence, relative to the current
 * state of the cDBG. Can either be: a null segment representing a portion
 * of the sequence already in the graph; a decision k-mer; or a new unitig or
 * portion of a unitig.
 */
struct compact_segment {
    // anchors are:
    //    1) if the segment is a unitig segment,
    //       the left or rightmost k-mer in the segment if it has no
    //       immediate existing unitig neighbor, or the unitig end it
    //       connects to if it does have that neighbor
    //    2) if the segment is a decision k-mer, the hash of that k-mer
    //       is both left and right anchor
    hash_t left_anchor;
    hash_t right_anchor;
    // has_left_unode and has_right_unode tell us if the immediate neighbors
    // of the segment are existung unitig ends at the time of creation. this
    // could be deduced using left_anchor and right_anchor along with the set
    // of new k-mers, but this saves a lot of lookups
    bool has_left_unode;
    bool has_right_unode;
    // whether this segment represents a decision k-mer
    bool is_decision_kmer;
    // start position of the segment within the originating sequence
    size_t start_pos;
    // length of the segment sequence (from beginning of first k-mer to
    // end of last k-mer)
    size_t length;

    // the default constructor creates a null segment
    compact_segment()
        : left_anchor(0),
          right_anchor(0),
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
        return (left_anchor == right_anchor) && !is_decision_kmer;
    }
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
    using AssemblerType::find_left_kmers;
    using AssemblerType::find_right_kmers;
    using AssemblerType::gather_left;
    using AssemblerType::gather_right;

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

        KmerIterator<AssemblerType> iter(sequence, this);
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

    bool get_decision_neighbors(const string& root,
                                NeighborBundle& result) {
        return get_decision_neighbors(this, root, result);
    }

    bool get_decision_neighbors(NeighborBundle& result) {
        return get_decision_neighbors(this, result);
    }

    bool get_decision_neighbors(AssemblerType* shifter,
                                const string& root,
                                NeighborBundle& result) {
        shifter->set_cursor(root);
        return get_decision_neighbors(shifter, result);
    }

    bool get_decision_neighbors(AssemblerType* shifter,
                                NeighborBundle& result) {

        auto left_kmers = shifter->find_left_kmers();
        auto right_kmers = shifter->find_right_kmers();
        
        if (left_kmers.size() > 1 || right_kmers.size() > 1) {
            result = std::make_pair(left_kmers, right_kmers);
            return true;
        } else {
            return false;
        }
    }

    compact_segment init_segment(hash_t left_anchor,
                                 bool has_left_unode,
                                 size_t start_pos) {
        compact_segment segment;
        segment.start_pos = start_pos;
        segment.has_left_unode = has_left_unode;
        segment.left_anchor = left_anchor;
        segment.is_decision_kmer = false;

        return segment;
    }

    void finish_segment(compact_segment& segment,
                        size_t end,
                        hash_t right_anchor,
                        bool has_right_unode,
                        deque<compact_segment>& segments) {
        
        segment.length = end - segment.start_pos + this->_K;
        segment.right_anchor = right_anchor;
        segment.has_right_unode = has_right_unode;
        segments.push_back(segment);
        pdebug("segment ended, start=" << segment.start_pos
               << " cur_pos=" << end);
    }

    void finish_decision_segment(compact_segment& segment,
                                deque<compact_segment>& segments) {
        segment.length = this->_K;
        segment.right_anchor = segment.left_anchor;
        segment.is_decision_kmer = true;
        segments.push_back(segment);
        segment.has_right_unode = false;
        segment.has_left_unode = false;

        pdebug("Created decision segment, start=" << segment.start_pos);       
    }

    void find_new_segments(const string& sequence,
                           set<hash_t>& new_kmers,
                           deque<compact_segment>& segments,
                           set<hash_t>& new_decision_kmers,
                           deque<NeighborBundle>& decision_neighbors
                           ) {


        vector<bool> kmer_new;
        vector<hash_t> hashes;
        dbg->add_sequence(sequence, hashes, kmer_new);
        pdebug("k-mers new: " << kmer_new);
        
        size_t pos = 0;
        hash_t prev_hash = hashes.front();
        bool cur_new = false, prev_new = false, is_decision = false;
        string kmer_seq;

        compact_segment current_segment;
        this->set_cursor(sequence);

        segments.push_back(compact_segment()); // place a null segment
        for (auto cur_hash : hashes) {
            cur_new = kmer_new[pos];

            if(cur_new) {
                new_kmers.insert(cur_hash);
                kmer_seq = sequence.substr(pos, this->_K);

                if(!prev_new || new_decision_kmers.count(prev_hash)) {
                    pdebug("old -> new, or prev d-kmer (pos=" << pos << ")");
                    this->set_cursor(kmer_seq);
                    hash_t left_anchor;
                    bool has_left_unode = false;
                    if (pos == 0) {
                        vector<shift_t> lneighbors = filter_nodes(gather_left());
                        if (lneighbors.size() == 1 &&
                            cdbg.has_unode_end(lneighbors.front().hash)) {
                            
                            left_anchor = lneighbors.front().hash;
                            has_left_unode = true;
                        } else {
                            left_anchor = cur_hash;
                        }
                    } else if (cdbg.has_unode_end(prev_hash)) {
                        left_anchor = prev_hash;
                        has_left_unode = true;
                    } else {
                        left_anchor = cur_hash;
                    }
                    current_segment = init_segment(left_anchor, has_left_unode, pos);

                } else {
                    this->shift_right(kmer_seq.back());
                }
                
                NeighborBundle neighbors;
                neighbors.first.clear();
                neighbors.second.clear();
                if ((is_decision = get_decision_neighbors(kmer_seq,
                                                          neighbors)) == true) {
                    decision_neighbors.push_back(neighbors);
                }
                if(is_decision) {
                    pdebug("new k-mer & decision" << this->get_cursor()
                           << ", " << kmer_seq);
                    new_decision_kmers.insert(cur_hash);

                    if (pos > 0 && prev_new && !new_decision_kmers.count(prev_hash)) {
                        finish_segment(current_segment,
                                       pos - 1,
                                       prev_hash,
                                       false,
                                       segments);
                    }

                    auto decision_segment = init_segment(cur_hash,
                                                         false,
                                                         pos);
                    finish_decision_segment(decision_segment,
                                            segments);
                }
            } else if (prev_new) {
                pdebug("new -> old");
                hash_t right_anchor;
                bool has_right_unode = false;
                if (cdbg.has_unode_end(cur_hash)) {
                    right_anchor = cur_hash;
                    has_right_unode = true;
                } else {
                    right_anchor = prev_hash;
                }
                finish_segment(current_segment, 
                               pos - 1,
                               right_anchor,
                               has_right_unode,
                               segments);
                segments.push_back(compact_segment()); // null segment
            } 

            ++pos;
            prev_hash = cur_hash;
            prev_new = cur_new;
        }

        // make sure we close current segment if necessary
        if (cur_new && !new_decision_kmers.count(hashes.back())) {
            hash_t right_anchor;
            bool has_right_unode = false;

            vector<shift_t> rneighbors = filter_nodes(gather_right());
            if (rneighbors.size() == 1 &&
                cdbg.has_unode_end(rneighbors.front().hash)) {
                right_anchor = rneighbors.front().hash;
                has_right_unode = true;
            } else {
                right_anchor = hashes.back();
            }
            pdebug("sequence ended on new k-mer");
            finish_segment(current_segment,
                           pos - 1, // we incr'd pos...
                           right_anchor,
                           has_right_unode,
                           segments);
        }

        if (cur_new) {
            segments.push_back(compact_segment()); // null segment
        }

        //pdebug("Segments: " << new_segment_sequences);
    }

    void update_from_segments(const string& sequence,
                              set<hash_t>& new_kmers,
                              deque<compact_segment>& segments,
                              set<hash_t> & new_decision_kmers,
                              deque<NeighborBundle>& decision_neighbors
                              ) {

        if (segments.size() == 0) {
            return;
        }

        // Pub: Algorithm 3: InsertSegment

        auto it = segments.begin();
        compact_segment& u = *it, v = *(++it), w = *(++it);

        while (it != segments.end()) {
            if (v.is_decision_kmer) {
                kmer_t decision_kmer(v.left_anchor,
                                     sequence.substr(v.start_pos, this->_K));
                _new_decision_kmer(decision_kmer,
                                   decision_neighbors.front(),
                                   new_kmers);
                _induce_decision_nodes(decision_kmer,
                                       decision_neighbors.front(),
                                       new_kmers);
                decision_neighbors.pop_front();
            } else if (!v.is_null()) {
                // v is a regular segment
                // if u is null, then we need to check for left induced d-nodes
                // or for a unitig connection
                if (u.is_null()) {
                    if (!v.has_left_unode) {
                        kmer_t root(v.left_anchor,
                                    sequence.substr(v.start_pos, this->_K));
                        // TODO: future optimization. we know root is not a decision
                        // k-mer, so in this case, unless v starts at pos = 0 in the sequence,
                        // the only neighbor is the preceeding k-mer in the sequence
                        _induce_decision_nodes_left(root, new_kmers);
                    }
                }

                // v is regular segment: if w is null, check for right induced
                // d-nodes or for a unitig connection
                if (w.is_null()) {
                    if (!v.has_right_unode) {
                        kmer_t root(v.right_anchor,
                                    sequence.substr(v.start_pos + v.length - this->_K, this->_K));
                        // TODO: future optimization. we know root is not a decision
                        // k-mer, so in this case, unless v starts at pos = 0 in the sequence,
                        // the only neighbor is the preceeding k-mer in the sequence
                        _induce_decision_nodes_right(root, new_kmers);
                    }
                }

                _update_unode(u, sequence);
            }

            u = v;
            v = w;
            w = *(++it);
        }
    }

    void update_sequence(const string& sequence) {
        set<hash_t> new_kmers;
        deque<compact_segment> segments;
        set<hash_t> new_decision_kmers;
        deque<NeighborBundle> decision_neighbors;

        find_new_segments(sequence,
                          new_kmers,
                          segments,
                          new_decision_kmers,
                          decision_neighbors);

        update_from_segments(sequence,
                             new_kmers,
                             segments,
                             new_decision_kmers,
                             decision_neighbors);
    }

    void _new_decision_kmer(kmer_t kmer,
                            NeighborBundle& neighbors,
                            set<hash_t>& neighbor_mask) {
        /* Handle a new k-mer that is also
         * a decision node.
         */
        cdbg.build_dnode_marker(kmer.hash);
        notify_build_dnode(kmer.hash, kmer.kmer);
    }

    void _induce_decision_nodes(kmer_t kmer,
                                NeighborBundle& neighbors,
                                set<hash_t>& neighbor_mask) {

        _induce_decision_nodes_left(kmer, neighbors, neighbor_mask);
        _induce_decision_nodes_right(kmer, neighbors, neighbor_mask);
    }

    uint8_t _induce_decision_nodes_left(kmer_t kmer,
                                        set<hash_t>& neighbor_mask) {

        this->set_cursor(kmer.kmer);
        NeighborBundle bundle;
        bundle.first = find_left_kmers();

        if (bundle.first.size()) {
            return _induce_decision_nodes_left(kmer,
                                               bundle,
                                               neighbor_mask);
        }
        return 0;
    }

    uint8_t _induce_decision_nodes_left(kmer_t kmer,
                                        NeighborBundle& neighbors,
                                        set<hash_t>& neighbor_mask) {

        // decision k-mers which are also new k-mers
        // cannot split existing unitigs. however, they can induce
        // an existing k-mer to be a decision k-mer, which can split
        // existing unitigs. so, we filter out neighbors of
        // the new decision k-mer which already exist and are already known
        // to the cDBG to be decision k-mers
                            

        vector<kmer_t> lneighbors;
        std::copy_if(neighbors.first.begin(),
                     neighbors.first.end(), 
                     std::back_inserter(lneighbors),
                     [&] (kmer_t neighbor) { return
                         !(neighbor_mask.count(neighbor.hash) ||
                          cdbg.query_dnode_marker(neighbor.hash));
                     });

        for (auto lneighbor : lneighbors) {
            NeighborBundle neighbors;
            if (get_decision_neighbors(lneighbor.kmer,
                                       neighbors)) {
                // induced decision k-mer

                pdebug("Induced d-node: " << lneighbor.hash << ", " << lneighbor.kmer);
                cdbg.build_dnode_marker(lneighbor.hash);
                notify_build_dnode(lneighbor.hash, lneighbor.kmer);
                _split_unode(lneighbor, neighbors, neighbor_mask); 
            }
            neighbors.first.clear();
            neighbors.second.clear();
        }

        return lneighbors.size();
    }

    uint8_t _induce_decision_nodes_right(kmer_t kmer,
                                         set<hash_t>& neighbor_mask) {

        this->set_cursor(kmer.kmer);
        NeighborBundle bundle;
        bundle.second = find_right_kmers();
        
        if (bundle.second.size()) {
            return _induce_decision_nodes_right(kmer,
                                                bundle,
                                                neighbor_mask);
        }
        return 0;
    }

    uint8_t _induce_decision_nodes_right(kmer_t kmer,
                                         NeighborBundle& neighbors,
                                         set<hash_t>& neighbor_mask) {

        // see _induce_decision_nodes_left for information

        vector<kmer_t> rneighbors;
        std::copy_if(neighbors.second.begin(),
                     neighbors.second.end(), 
                     std::back_inserter(rneighbors),
                     [&] (kmer_t neighbor) { return
                         !(neighbor_mask.count(neighbor.hash) ||
                          cdbg.query_dnode_marker(neighbor.hash));
                     });

        for (auto rneighbor : rneighbors) {
            NeighborBundle neighbors;
            if (get_decision_neighbors(rneighbor.kmer,
                                       neighbors)) {
                // induced decision k-mer
                pdebug("Induced d-node: " << rneighbor.hash << ", " << rneighbor.kmer);
                cdbg.build_dnode_marker(rneighbor.hash);
                notify_build_dnode(rneighbor.hash, rneighbor.kmer);
                _split_unode(rneighbor, neighbors, neighbor_mask); 
            }
            neighbors.first.clear();
            neighbors.second.clear();
        }

        return rneighbors.size();
    }

    void _split_unode(kmer_t root,
                      NeighborBundle& neighbors,
                      set<hash_t>& mask) {

        
    }

    void _update_unode(compact_segment& segment,
                       const string& sequence) {

    }
};

}
#undef pdebug
#endif
