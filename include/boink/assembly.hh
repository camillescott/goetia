/* assembly.hh -- boink traversal and assembly
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_ASSEMBLY_HH
#define BOINK_ASSEMBLY_HH

#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/kmeriterator.hh"

#include <deque>
#include <memory>
#include <set>
#include <unordered_set>
#include <vector>


# ifdef DEBUG_ASMLY
#   define pdebug(x) do { std::cerr << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif


namespace boink {

using boink::hashing::hash_t;
using boink::hashing::kmer_t;
using boink::hashing::shift_t;


typedef std::deque<char> Path;
typedef std::vector<std::string> StringVector;
typedef std::vector<kmer_t> KmerVector;

typedef std::pair<KmerVector, KmerVector> NeighborBundle;
typedef std::pair<kmer_t, NeighborBundle> DecisionKmer;


struct DecisionKmerHash {
    // Just uses the hash of the kmer_t from the pair
    inline hash_t operator()(const DecisionKmer& n) const {
        return n.first.hash;
    }
};

struct DecisionKmerComp {
    bool operator()(const DecisionKmer& a, const DecisionKmer& b) {
        return a.first.hash < b.first.hash;
    }
};

typedef std::set<DecisionKmer, DecisionKmerComp> DecisionKmerSet;
typedef std::unordered_set<DecisionKmer, DecisionKmerHash> DecisionKmerUSet;

namespace TraversalState {

    enum State {
        STOP_FWD,
        STOP_RC,
        DECISION_FWD,
        DECISION_RC,
        STOP_SEEN,
        STOP_MASKED,
        BAD_SEED
    };
}

template <class GraphType>
struct Traverse {

    using BaseShifter = typename GraphType::shifter_type;
    
    typedef GraphType graph_type;
    typedef BaseShifter shifter_type;


    typedef TraversalState::State State;
    typedef std::pair<State, hashing::hash_t> EndState;

    class dBG: public BaseShifter {

    protected:

        std::set<hash_t> seen;

    public:

        using BaseShifter::set_cursor;
        using BaseShifter::get_cursor;
        using BaseShifter::get;
        using BaseShifter::hash;
        using BaseShifter::gather_left;
        using BaseShifter::gather_right;
        using BaseShifter::shift_left;
        using BaseShifter::shift_right;

        std::shared_ptr<GraphType> graph;

        dBG(BaseShifter const& shifter)
            : BaseShifter(shifter)
        {
        }

        dBG(uint16_t K)
            : BaseShifter(K)
        {
        }

        static std::shared_ptr<dBG> build(uint16_t K) {
            return std::make_shared<dBG>(K);
        }

        static std::shared_ptr<dBG> build(BaseShifter const& shifter) {
            return std::make_shared<dBG>(shifter);
        }

        void clear_seen() {
            seen.clear();
        }

        size_t degree_left(GraphType * graph) {
            auto neighbors = this->gather_left();
            return count_nodes(graph, neighbors);
        }

        size_t degree_right(GraphType * graph) {
            auto neighbors = this->gather_right();
            return count_nodes(graph, neighbors);
        }

        size_t degree(GraphType * graph) {
            return degree_left(graph) + degree_right(graph);
        }

        size_t degree_left(GraphType * graph,
                           std::set<hash_t>& extras) {
            auto neighbors = this->gather_left();
            return count_nodes(graph, neighbors, extras);
        }

        size_t degree_right(GraphType * graph,
                            std::set<hash_t>& extras) {
            auto neighbors = this->gather_right();
            return count_nodes(graph, neighbors, extras);
        }

        size_t degree(GraphType * graph,
                      std::set<hash_t>& extras) {
            return degree_left(graph, extras) + degree_right(graph, extras);
        }

        size_t get_left(GraphType * graph,
                        shift_t& result) {
            std::vector<shift_t> neighbors = this->gather_left();
            auto n_left = reduce_nodes(graph, neighbors, result);

            if (n_left == 1) {
                this->shift_left(result.symbol);
            }
            return n_left;
        }

        size_t get_right(GraphType * graph,
                         shift_t& result) {
            std::vector<shift_t> neighbors = this->gather_right();
            auto n_right = reduce_nodes(graph, neighbors, result);

            if (n_right == 1) {
                this->shift_right(result.symbol);
            }
            return n_right;
        }


        size_t count_nodes(GraphType *                 graph,
                           const std::vector<shift_t>& nodes) {
            uint8_t n_found = 0;
            for (auto node: nodes) {
                if(graph->query(node.hash)) {
                    ++n_found;
                }
            }
            return n_found;
        }

        size_t count_nodes(GraphType *                 graph,
                           const std::vector<shift_t>& nodes,
                           std::set<hash_t>&           extras) {
            uint8_t n_found = 0;
            for (auto node: nodes) {
                if(graph->query(node.hash) ||
                   extras.count(node.hash)) {
                    ++n_found;
                }
            }
            return n_found;
        }

        size_t reduce_nodes(GraphType *                 graph,
                            const std::vector<shift_t>& nodes,
                            shift_t&                    result) {
            uint8_t n_found = 0;
            for (auto node : nodes) {
                if(graph->query(node.hash)) {
                    ++n_found;
                    if (n_found > 1) {
                        return n_found;
                    }
                    result = node;
                }
            }
            return n_found;
        }

        size_t reduce_nodes(GraphType *                 graph,
                            const std::vector<shift_t>& nodes,
                             shift_t&                   result,
                             std::set<hash_t>&          extra) {
            uint8_t n_found = 0;
            for (auto node : nodes) {
                if(graph->query(node.hash) ||
                   extra.count(node.hash)) {
                    ++n_found;
                    if (n_found > 1) {
                        return n_found;
                    }
                    result = node;
                }
            }
            return n_found;
        }

        std::vector<shift_t> filter_nodes(GraphType * graph,
                                          const std::vector<shift_t>& nodes) {
            std::vector<shift_t> result;
            for (auto node : nodes) {
                if (graph->query(node.hash)) {
                    result.push_back(node);
                }
            }
            return result;
        }

        std::vector<shift_t> filter_nodes(GraphType *                 graph,
                                          const std::vector<shift_t>& nodes,
                                          std::set<hash_t>&           extra) {
            std::vector<shift_t> result;
            for (auto node : nodes) {
                if (graph->query(node.hash) ||
                    extra.count(node.hash)) {
                    result.push_back(node);
                }
            }
            return result;
        }

        std::vector<kmer_t> find_left_kmers(GraphType * graph) {
            auto root = this->get_cursor();
            auto filtered = filter_nodes(graph, this->gather_left());
            return graph->build_left_kmers(filtered, root);
        }

        std::vector<kmer_t> find_right_kmers(GraphType * graph) {
            auto root = this->get_cursor();
            auto filtered = filter_nodes(graph, this->gather_right());
            return graph->build_right_kmers(filtered, root);
        }

        std::vector<kmer_t> find_left_kmers(GraphType *       graph,
                                            std::set<hash_t>& extras) {
            auto root = this->get_cursor();
            auto filtered = filter_nodes(graph, this->gather_left(), extras);
            return graph->build_left_kmers(filtered, root);
        }

        std::vector<kmer_t> find_right_kmers(GraphType *       graph,
                                             std::set<hash_t>& extras) {
            auto root = this->get_cursor();
            auto filtered = filter_nodes(graph, this->gather_right(), extras);
            return graph->build_right_kmers(filtered, root);
        }

        EndState traverse_left(GraphType *        graph,
                               const std::string& seed,
                               Path&              path,
                               std::set<hash_t>&  mask) {

            this->set_cursor(seed);
            if (!graph->query(this->get())) {
                return {State::BAD_SEED, this->get()};
            } 
            this->get_cursor(path);
            return traverse_left(graph, path, mask);
        }

        EndState traverse_left(GraphType * graph,
                               Path&       path,
                               std::set<hash_t>& mask) {

            auto end_hash = this->get();
            this->seen.clear();
            this->seen.insert(this->get());

            shift_t next;
            uint8_t n_left;
            while (1) {
                if (degree_right(graph) > 1) {
                    pdebug("Stop: reverse d-node");
                    path.pop_front();
                    return {State::DECISION_RC, end_hash};
                }

                n_left = this->reduce_nodes(graph,
                                            this->gather_left(),
                                            next);
                end_hash = this->get();

                if (n_left > 1) {
                    pdebug("Stop: forward d-node");
                    return {State::DECISION_FWD, end_hash};
                }

                if (n_left == 0) {
                    pdebug("Stop: no neighbors.");
                    return {State::STOP_FWD, end_hash};
                }

                if (this->seen.count(next.hash)) {
                    pdebug("Stop: hit seen k-mer.");
                    return {State::STOP_SEEN, end_hash};       
                       
                }
                if (mask.count(next.hash)) {
                    pdebug("Stop: hit masked k-mer.");
                    return {State::STOP_MASKED, end_hash};     
                }
                
                this->shift_left(next.symbol);
                path.push_front(next.symbol);
                this->seen.insert(next.hash);
            }
        }

        EndState traverse_right(GraphType *        graph,
                                const std::string& seed,
                                Path&              path,
                                std::set<hash_t>&  mask) {

            this->set_cursor(seed);
            if (!graph->query(this->get())) {
                return {State::BAD_SEED, this->get()};
            }
            this->get_cursor(path);
            return traverse_right(graph, path, mask);
        }

        EndState traverse_right(GraphType *       graph,
                                Path&             path,
                                std::set<hash_t>& mask) {

            auto end_hash = this->get();
            this->seen.clear();
            this->seen.insert(this->get());
            
            shift_t next;
            uint8_t n_right;
            while (1) {
                if (degree_left(graph) > 1) {
                    path.pop_back();
                    return {State::DECISION_RC, end_hash};
                }

                n_right = this->reduce_nodes(graph,
                                             this->gather_right(),
                                             next);
                end_hash = this->get();

                if (n_right > 1) {
                    return {State::DECISION_FWD, end_hash};
                }

                if (n_right == 0) {
                    return {State::STOP_FWD, end_hash};
                }

                if (this->seen.count(next.hash)) {
                    pdebug("Stop: hit seen k-mer.");
                    return {State::STOP_SEEN, end_hash};       
                       
                }
                if (mask.count(next.hash)) {
                    pdebug("Stop: hit masked k-mer.");
                    return {State::STOP_MASKED, end_hash};     
                }

                this->shift_right(next.symbol);
                path.push_back(next.symbol);
                this->seen.insert(next.hash);
            }
        }

        std::pair<EndState, EndState> traverse(GraphType *        graph,
                                               const std::string& seed,
                                               Path&              path,
                                               std::set<hash_t>   mask) {

            this->set_cursor(seed);
            if (!graph->query(this->hash(seed))) {
                return {{State::BAD_SEED, this->get()},
                        {State::BAD_SEED, this->get()}};
            }
            this->get_cursor(path);
            auto state_left = traverse_left(graph, path, mask);
            this->set_cursor(seed);
            auto state_right = traverse_right(graph, path, mask);
            return {state_left, state_right};
        }
                      
        std::string to_string(Path& path) {
            return std::string(path.begin(), path.end());
        }

        bool is_decision_kmer(GraphType *        graph,
                              const std::string& node,
                              uint8_t&           degree) {
            this->set_cursor(node);
            return is_decision_kmer(graph, degree);
        }

        bool is_decision_kmer(GraphType *        graph,
                              const std::string& node) {
            this->set_cursor(node);
            return this->degree_left(graph) > 1 || this->degree_right(graph) > 1;
        }

        bool is_decision_kmer(GraphType * graph,
                              uint8_t&    degree) {
            uint8_t ldegree, rdegree;
            ldegree = this->degree_left(graph);
            rdegree = this->degree_right(graph);
            degree = ldegree + rdegree;
            return ldegree > 1 || rdegree > 1;
        }

        void find_decision_kmers(GraphType *                  graph,
                                 const std::string&           sequence,
                                 std::vector<uint32_t>&       decision_positions,
                                 std::vector<hash_t>&         decision_hashes,
                                 std::vector<NeighborBundle>& decision_neighbors) {

            hashing::KmerIterator<dBG> iter(sequence, this);
            size_t pos = 0;
            while(!iter.done()) {
                hash_t h = iter.next();
                NeighborBundle neighbors;
                if (get_decision_neighbors(graph,
                                           iter.shifter,
                                           neighbors)) {

                    decision_neighbors.push_back(neighbors);
                    decision_positions.push_back(pos);
                    decision_hashes.push_back(h);
                }
            
                ++pos;
           }
        }

        bool get_decision_neighbors(GraphType *        graph,
                                    const std::string& root,
                                    NeighborBundle&    result,
                                    std::set<hash_t>&  union_nodes) {

            return get_decision_neighbors(graph, this, root, result, union_nodes);
        }

        bool get_decision_neighbors(GraphType *       graph,
                                    NeighborBundle&   result,
                                    std::set<hash_t>& union_nodes) {
            return get_decision_neighbors(graph, this, result, union_nodes);
        }

        bool get_decision_neighbors(GraphType *        graph,
                                    dBG *              shifter,
                                    const std::string& root,
                                    NeighborBundle&    result,
                                    std::set<hash_t>&  union_nodes) {
            shifter->set_cursor(root);
            return get_decision_neighbors(graph, shifter, result, union_nodes);
        }

        bool get_decision_neighbors(GraphType *       graph,
                                    dBG *             shifter,
                                    NeighborBundle&   result,
                                    std::set<hash_t>& union_nodes) {

            auto left_kmers = shifter->find_left_kmers(graph, union_nodes);
            auto right_kmers = shifter->find_right_kmers(graph, union_nodes);
            
            if (left_kmers.size() > 1 || right_kmers.size() > 1) {
                result = std::make_pair(left_kmers, right_kmers);
                return true;
            } else {
                return false;
            }
        }

        bool get_decision_neighbors(GraphType *        graph,
                                    const std::string& root,
                                    NeighborBundle&    result) {
            return get_decision_neighbors(graph, this, root, result);
        }

        bool get_decision_neighbors(GraphType *     graph,
                                    NeighborBundle& result) {
            return get_decision_neighbors(graph, this, result);
        }

        bool get_decision_neighbors(GraphType *        graph,
                                    dBG *              shifter,
                                    const std::string& root,
                                    NeighborBundle&    result) {
            shifter->set_cursor(root);
            return get_decision_neighbors(graph, shifter, result);
        }

        bool get_decision_neighbors(GraphType *     graph,
                                    dBG *           shifter,
                                    NeighborBundle& result) {

            auto left_kmers = shifter->find_left_kmers(graph);
            auto right_kmers = shifter->find_right_kmers(graph);
            
            if (left_kmers.size() > 1 || right_kmers.size() > 1) {
                result = std::make_pair(left_kmers, right_kmers);
                return true;
            } else {
                return false;
            }
        }

    };
};
}

#undef pdebug
#endif
