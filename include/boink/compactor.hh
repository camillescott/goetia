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

namespace boink {

#include "boink/assembly.hh"
#include "boink/hashing.hh"
#include "boink/dbg.hh"
#include "boink/cdbg.hh"
#include "boink/minimizers.hh"


# ifdef DEBUG_COMPACTOR
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


template <class GraphType>
class StreamingCompactor : public AssemblerMixin<GraphType> {

protected:

    uint64_t _minimizer_window_size;

public:

    GraphType * dbg;
    cDBG cdbg;

    using AssemblerType = AssemblerMixin<GraphType>;
    using AssemblerType::seen;
    using AssemblerType::get_left;
    using AssemblerType::get_right;
    using AssemblerType::degree_left;
    using AssemblerType::degree_right;
    using AssemblerType::count_nodes;
    using AssemblerType::filter_nodes;

    using ShifterType = typename GraphType::shifter_type;

    StreamingCompactor(GraphType * dbg,
                       uint64_t minimizer_window_size=8) :
        AssemblerMixin<GraphType>(dbg),
        _minimizer_window_size(minimizer_window_size),
        dbg(dbg),
        cdbg(dbg->K()) {

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

    bool insert_sequence(const string& sequence,
                         vector<uint32_t>& decision_positions,
                         HashVector& decision_hashes,
                         vector<NeighborBundle>& decision_neighbors) {

        if (!dbg->add_sequence(sequence)) {
            return false; // if there were no new k-mers, nothing to do
        } else {
            find_decision_kmers(sequence,
                                decision_positions,
                                decision_hashes,
                                decision_neighbors);
            return true;
        }
    }

    void update_sequence(const string& sequence) {
        if(!dbg->add_sequence(sequence)) {
            return;
        } else {
            update_cdbg(sequence);
        }
    }

    void update_cdbg(const string& sequence) {
        vector<DecisionNode*> disturbed_dnodes;
        vector<NeighborBundle> disturbed_neighbors;
        find_disturbed_dnodes(sequence,
                              disturbed_dnodes,
                              disturbed_neighbors);
        if (disturbed_dnodes.size() == 0) {
            _update_linear(sequence);
        } else {
            _update_from_dnodes(sequence,
                                disturbed_dnodes,
                                disturbed_neighbors);
        }
    }

    void _update_linear(const string& sequence) {
        pdebug("Update linear");
        this->set_cursor(sequence.substr(0, this->_K));
        Path segment;
        segment.insert(segment.end(), sequence.begin(), sequence.end());
        this->compactify_left(segment);
        hash_t stopped_at = this->get();
        char after = *(segment.begin()+this->_K);
        junction_t left_junc = make_pair(stopped_at,
                                         this->shift_right(after));


        this->set_cursor(sequence.substr(sequence.length()-this->_K,
                                         this->_K));
        this->compactify_right(segment);
        stopped_at = this->get();
        after = *(segment.end()-(this->_K)-1);
        junction_t right_junc = make_pair(this->shift_left(after),
                               stopped_at);

        string segment_seq = this->to_string(segment);
        HashVector tags;
        WKMinimizer<ShifterType> minimizer(_minimizer_window_size,
                                           this->_K);
        tags = minimizer.get_minimizer_values(segment_seq);
        cdbg.build_unode(tags, segment_seq, left_junc, right_junc);
    }

    void _update_from_dnodes(const string& sequence,
                             vector<DecisionNode*>& disturbed_dnodes,
                             vector<NeighborBundle>& disturbed_neighbors) {

        std::set<hash_t> updated_dnodes;
        std::set<junction_t> updated_junctions;
        pdebug(disturbed_dnodes.size() << " on d-node queue, " <<
               disturbed_neighbors.size() << " on neighbor queue");

        while(!disturbed_dnodes.empty()) {
            DecisionNode * root_dnode = disturbed_dnodes.back();
            NeighborBundle root_neighbors = disturbed_neighbors.back();
            disturbed_dnodes.pop_back();
            disturbed_neighbors.pop_back();

            pdebug("updating from " << *root_dnode);

            if (updated_dnodes.count(root_dnode->node_id)) {
                continue;
            }

            for (kmer_t left_neighbor : root_neighbors.first) {
                junction_t right_junc = make_pair(left_neighbor.hash,
                                                  root_dnode->node_id);
				pdebug("left neighbor: " << left_neighbor.kmer 
                        << " junction: " << right_junc);
                if (updated_junctions.count(right_junc)) {
                    pdebug("Already updated from, continuing");
                    continue;
                }
                
                HashVector tags;
                std::string segment_seq;
                junction_t left_junc;
                if (cdbg.get_dnode(left_neighbor.hash) != nullptr) {
                    // trivial unode
                    left_junc = right_junc;
                    segment_seq = left_neighbor.kmer
                                  + root_dnode->sequence.back();
                    pdebug("Trivial u-node, junctions " << left_junc
                           << ", " << right_junc << ", sequence="
                           << segment_seq);
                } else {
                    this->set_cursor(left_neighbor.kmer);
                    Path this_segment;
                    this->get_cursor(this_segment);
                    pdebug("Start assembly left at " << this->get_cursor());
                    this_segment.push_back(root_dnode->sequence.back());

                    InteriorMinimizer<hash_t> minimizer(_minimizer_window_size);
                    this->compactify_left(this_segment, minimizer);
                    
                    hash_t stopped_at = this->get();
                    std::string after = std::string(this_segment.begin()+1,
                                                    this_segment.begin()+this->_K+1);
                    left_junc = make_pair(this->get(),
                                          this->shift_right(after.back()));
                    segment_seq = this->to_string(this_segment);
                    pdebug("Assembled left, stopped left of " << after
                            << ", shifting right on " << after.back()
                            << " with junction " << left_junc
                            << ", sequence=" << segment_seq);
                    tags = minimizer.get_minimizer_values();
                }
                updated_junctions.insert(left_junc);
                updated_junctions.insert(right_junc);

                cdbg.build_unode(tags, segment_seq, left_junc, right_junc);
            }

            for (kmer_t right_neighbor : root_neighbors.second) {
                junction_t left_junc = make_pair(root_dnode->node_id,
                                                 right_neighbor.hash);
				pdebug("right neighbor: " << right_neighbor.kmer 
                        << " junction: " << left_junc);
                if (updated_junctions.count(left_junc)) {
                    pdebug("Already updated from, continuing");
                    continue;
                }

                HashVector tags;
                std::string segment_seq;
                junction_t right_junc;
                if (cdbg.get_dnode(right_neighbor.hash) != nullptr) {
                    // trivial unode
                    right_junc = left_junc;
                    segment_seq = root_dnode->sequence.front()
                                  + right_neighbor.kmer;
                    pdebug("Trivial u-node, junctions " << left_junc
                           << ", " << right_junc << ", sequence="
                           << segment_seq);
                } else {
                    this->set_cursor(right_neighbor.kmer);
                    Path this_segment;
                    this->get_cursor(this_segment);
                    this_segment.push_front(root_dnode->sequence.front());

                    InteriorMinimizer<hash_t> minimizer(_minimizer_window_size);
                    this->compactify_right(this_segment, minimizer);

                    hash_t stopped_at = this->get();
                    std::string after = std::string(this_segment.end()-(this->_K)-1,
                                                    this_segment.end()-1);
                    right_junc = make_pair(this->shift_left(after.front()),
                                                      stopped_at);
                    segment_seq = this->to_string(this_segment);
                    pdebug("Assembled right, stopped right of " << after
                            << " with junction " << right_junc
                            << ", sequence=" << segment_seq);
                    tags = minimizer.get_minimizer_values();
                }
                updated_junctions.insert(left_junc);
                updated_junctions.insert(right_junc);

                cdbg.build_unode(tags, segment_seq, left_junc, right_junc);
            }

            updated_dnodes.insert(root_dnode->node_id);
        }
        pdebug("FINISHED _update_from_dnodes with " << updated_dnodes.size()
                << " updated from; " << cdbg.n_decision_nodes() 
                << " total d-nodes  and " << cdbg.n_unitig_nodes() << " u-nodes");
    }

    bool is_decision_kmer(uint8_t& degree) {
        // oops better put notes here
        uint8_t ldegree, rdegree;
        ldegree = this->degree_left();
        rdegree = this->degree_right();
        degree = ldegree + rdegree;
        return ldegree > 1 || rdegree > 1;
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

    void find_disturbed_dnodes(const string& sequence,
                               vector<DecisionNode*>& disturbed_dnodes,
                               vector<NeighborBundle>& disturbed_neighbors) {
        
        KmerIterator<typename GraphType::shifter_type> iter(sequence, this->_K);

        // first we have to search for induced decision nodes in-incident to
        // the first k-mer
        vector<shift_t> left_neighbors = filter_nodes(iter.shifter.gather_left());
        typename GraphType::shifter_type flank_shifter(iter.shifter);
        const char suffix = sequence[this->_K-1];
        for (shift_t left_neighbor : left_neighbors) {
            flank_shifter.shift_left(left_neighbor.symbol);
            NeighborBundle decision_neighbors;
            string kmer = flank_shifter.get_cursor();

            if (get_decision_neighbors(flank_shifter,
                                       kmer,
                                       decision_neighbors)) {
                DecisionNode * dnode;
                dnode = cdbg.build_dnode(flank_shifter.get(), kmer);
                disturbed_dnodes.push_back(dnode);
                disturbed_neighbors.push_back(decision_neighbors);
            }
            flank_shifter.shift_right(suffix);
        }

        size_t pos = 0;
        while(!iter.done()) {
            hash_t h = iter.next();
            NeighborBundle decision_neighbors;
            string kmer = sequence.substr(pos, this->_K);
            if (get_decision_neighbors(iter.shifter,
                                       kmer,
                                       decision_neighbors)) {
                pdebug("Found d-node " << h << ", " << kmer <<
                       " ldegree " << decision_neighbors.first.size() <<
                       " rdegree " << decision_neighbors.second.size());
                DecisionNode * dnode;
                dnode = cdbg.build_dnode(h, kmer);
                disturbed_dnodes.push_back(dnode);
                disturbed_neighbors.push_back(decision_neighbors);
            }
            ++pos;
        }

        // and get the right flanking nodes
        vector<shift_t> right_neighbors = filter_nodes(iter.shifter.gather_right());
        const char prefix = sequence[sequence.length() - this->_K];
        for (shift_t right_neighbor : right_neighbors) {
            iter.shifter.shift_right(right_neighbor.symbol);
            NeighborBundle decision_neighbors;
            string kmer = iter.shifter.get_cursor();

            if (get_decision_neighbors(iter.shifter,
                                       kmer,
                                       decision_neighbors)) {
                DecisionNode * dnode;
                dnode = cdbg.build_dnode(iter.shifter.get(), kmer);
                disturbed_dnodes.push_back(dnode);
                disturbed_neighbors.push_back(decision_neighbors);
            }
            iter.shifter.shift_left(prefix);
        }
    }


};




}
#undef pdebug
#endif
