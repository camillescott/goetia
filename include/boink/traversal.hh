/**
 * (c) Camille Scott, 2019
 * File   : assembly.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 23.07.2019
 */


#ifndef BOINK_ASSEMBLY_HH
#define BOINK_ASSEMBLY_HH

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

typedef std::deque<char> Path;
typedef std::vector<std::string> StringVector;

/*
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
*/

namespace TraversalState {

    enum State {
        STOP_FWD,
        STOP_RC,
        DECISION_FWD,
        DECISION_RC,
        STOP_SEEN,
        STOP_MASKED,
        BAD_SEED,
        GRAPH_ERROR,
        STEP
    };
}


template <class GraphType>
struct Traverse {

    using BaseShifter = typename GraphType::shifter_type;
    
    typedef GraphType                         graph_type;
    typedef BaseShifter                       shifter_type;

    typedef typename shifter_type::hash_type  hash_type;
    typedef typename shifter_type::shift_type shift_type;
    typedef typename shifter_type::kmer_type  kmer_type;

    typedef std::pair<std::vector<kmer_type>,
                      std::vector<kmer_type>> NeighborBundle;

    typedef TraversalState::State       State;
    typedef std::pair<State, hash_type> EndState;

    class dBG: public BaseShifter {

    protected:

        std::set<hash_type> seen;

    public:

        using BaseShifter::set_cursor;
        using BaseShifter::get_cursor;
        using BaseShifter::get;
        using BaseShifter::hash;
        using BaseShifter::gather_left;
        using BaseShifter::gather_right;
        using BaseShifter::shift_left;
        using BaseShifter::shift_right;

        template<typename Arg1,
                 typename... Args>
        explicit dBG(Arg1 arg1, Args&&... args)
            : BaseShifter(arg1, std::forward<Args>(args)...)
        {
        }

        explicit dBG(BaseShifter const& shifter)
            : BaseShifter(shifter)
        {
        }

        template<typename... Args>
        static std::shared_ptr<dBG> build(Args&&... args) {
            return std::make_shared<dBG>(std::forward<Args>(args)...);
        }

        static std::shared_ptr<dBG> build(BaseShifter const& shifter) {
            return std::make_shared<dBG>(shifter);
        }

        void clear_seen() {
            seen.clear();
        }

        size_t in_degree(GraphType * graph) {
            auto neighbors = this->gather_left();
            return count_nodes(graph, neighbors);
        }

        size_t out_degree(GraphType * graph) {
            auto neighbors = this->gather_right();
            return count_nodes(graph, neighbors);
        }

        size_t degree(GraphType * graph) {
            return in_degree(graph) + out_degree(graph);
        }

        size_t in_degree(GraphType *          graph,
                           std::set<hash_type>& extras) {

            auto neighbors = this->gather_left();
            return count_nodes(graph, neighbors, extras);
        }

        size_t out_degree(GraphType *          graph,
                            std::set<hash_type>& extras) {

            auto neighbors = this->gather_right();
            return count_nodes(graph, neighbors, extras);
        }

        size_t degree(GraphType *          graph,
                      std::set<hash_type>& extras) {
            return in_degree(graph, extras) + out_degree(graph, extras);
        }

        size_t count_nodes(GraphType *                    graph,
                           const std::vector<shift_type>& nodes);

        size_t count_nodes(GraphType *                    graph,
                           const std::vector<shift_type>& nodes,
                           std::set<hash_type>&           extras);

        size_t reduce_nodes(GraphType *                    graph,
                            const std::vector<shift_type>& nodes,
                            shift_type&                    result);

        size_t reduce_nodes(GraphType *                    graph,
                            const std::vector<shift_type>& nodes,
                            shift_type&                    result,
                            std::set<hash_type>&           extra);

        std::vector<shift_type> filter_nodes(GraphType *                    graph,
                                             const std::vector<shift_type>& nodes);

        std::vector<shift_type> filter_nodes(GraphType *                    graph,
                                             const std::vector<shift_type>& nodes,
                                             std::set<hash_type>&           extra);

        std::pair<std::vector<shift_type>,
                  std::vector<shift_type>> filter_nodes(GraphType * graph,
                                                        const std::pair<std::vector<shift_type>,
                                                                        std::vector<shift_type>>& nodes);

        std::pair<std::vector<shift_type>,
                  std::vector<shift_type>> filter_nodes(GraphType *                               graph,
                                                        const std::pair<std::vector<shift_type>,
                                                                        std::vector<shift_type>>& nodes,
                                                        std::set<hash_type>&                      extras);

        std::vector<shift_type> in_neighbors(GraphType * graph) {
            auto root = this->get_cursor();
            return filter_nodes(graph, this->gather_left());
        }

        std::vector<shift_type> out_neighbors(GraphType * graph) {
            auto root = this->get_cursor();
            return filter_nodes(graph, this->gather_right());
        }

        std::pair<std::vector<shift_type>,
                  std::vector<shift_type>> neighbors(GraphType * graph) {
            auto _in = in_neighbors(graph);
            auto _out = out_neighbors(graph);
            return std::make_pair(_in, _out);
        }

        std::vector<kmer_type> find_left_kmers(GraphType * graph) {
            auto root = this->get_cursor();
            auto filtered = filter_nodes(graph, this->gather_left());
            return graph->build_left_kmers(filtered, root);
        }

        std::vector<kmer_type> find_right_kmers(GraphType * graph) {
            auto root = this->get_cursor();
            auto filtered = filter_nodes(graph, this->gather_right());
            return graph->build_right_kmers(filtered, root);
        }

        std::vector<kmer_type> find_left_kmers(GraphType *       graph,
                                               std::set<hash_type>& extras) {
            auto root = this->get_cursor();
            auto filtered = filter_nodes(graph, this->gather_left(), extras);
            return graph->build_left_kmers(filtered, root);
        }

        std::vector<kmer_type> find_right_kmers(GraphType *       graph,
                                                std::set<hash_type>& extras) {
            auto root = this->get_cursor();
            auto filtered = filter_nodes(graph, this->gather_right(), extras);
            return graph->build_right_kmers(filtered, root);
        }

        State look_state(std::vector<shift_type>& neighbors) {
            if (neighbors.size() > 1) {
                pdebug("Stop: forward d-node");
                return State::DECISION_FWD;
            } else if (neighbors.size() == 0) {
                pdebug("Stop: no neighbors.");
                return State::STOP_FWD;
            } else if (this->seen.count(neighbors.front().hash)) {
                pdebug("Stop: hit seen k-mer.");
                return State::STOP_SEEN;
            } else {
                return State::STEP;
            }
        }

        std::pair<State, std::vector<shift_type>> step_left(GraphType * graph) {
            auto neighbors = this->filter_nodes(graph, this->gather_left());
            auto state = look_state(neighbors);
            if (state == State::STEP) {
                this->shift_left(neighbors.front().symbol);
                this->seen.insert(neighbors.front().hash);
            }
            return {state, std::move(neighbors)};
        }

        std::pair<State, std::vector<shift_type>> step_left(GraphType * graph,
                                                            std::set<hash_type>& mask) {

            auto neighbors = this->filter_nodes(graph, this->gather_left());
            auto state = look_state(neighbors);
            if (state == State::STEP) {
                if (mask.count(neighbors.front().hash)) {
                    state = State::STOP_MASKED;
                } else {
                    this->shift_left(neighbors.front().symbol);
                    this->seen.insert(neighbors.front().hash);
                }
            }
            return {state, std::move(neighbors)};
        }

        std::pair<State, std::vector<shift_type>> step_right(GraphType * graph) {

            auto neighbors = this->filter_nodes(graph, this->gather_right());
            auto state = look_state(neighbors);
            if (state == State::STEP) {
                this->shift_right(neighbors.front().symbol);
                this->seen.insert(neighbors.front().hash);
            }
            return {state, std::move(neighbors)};
        }

        std::pair<State, std::vector<shift_type>> step_right(GraphType * graph,
                                                             std::set<hash_type>& mask) {
            auto neighbors = this->filter_nodes(graph, this->gather_right());
            auto state = look_state(neighbors);
            if (state == State::STEP) {
                if (mask.count(neighbors.front().hash)) {
                    state = State::STOP_MASKED;
                } else {
                    this->shift_right(neighbors.front().symbol);
                    this->seen.insert(neighbors.front().hash);
                }
            }
            return {state, std::move(neighbors)};
        }

        EndState walk_left(GraphType *          graph,
                               const std::string&   seed,
                               Path&                path,
                               std::set<hash_type>& mask);

        EndState walk_left(GraphType *          graph,
                               Path&                path,
                               std::set<hash_type>& mask);

        EndState walk_right(GraphType *          graph,
                                const std::string&   seed,
                                Path&                path,
                                std::set<hash_type>& mask);

        EndState walk_right(GraphType *       graph,
                                Path&             path,
                                std::set<hash_type>& mask);

        std::pair<EndState, EndState> walk(GraphType *         graph,
                                               const std::string&  seed,
                                               Path&               path,
                                               std::set<hash_type> mask);
                      
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
            return this->in_degree(graph) > 1 || this->out_degree(graph) > 1;
        }

        bool is_decision_kmer(GraphType * graph,
                              uint8_t&    degree) {

            uint8_t ldegree, rdegree;
            ldegree = this->in_degree(graph);
            rdegree = this->out_degree(graph);
            degree = ldegree + rdegree;
            return ldegree > 1 || rdegree > 1;
        }

        void find_decision_kmers(GraphType *                  graph,
                                 const std::string&           sequence,
                                 std::vector<uint32_t>&       decision_positions,
                                 std::vector<hash_type>&      decision_hashes,
                                 std::vector<NeighborBundle>& decision_neighbors) {

            hashing::KmerIterator<dBG> iter(sequence, this);
            size_t pos = 0;
            while(!iter.done()) {
                hash_type h = iter.next();
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

        bool get_decision_neighbors(GraphType *          graph,
                                    const std::string&   root,
                                    NeighborBundle&      result,
                                    std::set<hash_type>& union_nodes) {

            return get_decision_neighbors(graph, this, root, result, union_nodes);
        }

        bool get_decision_neighbors(GraphType *          graph,
                                    NeighborBundle&      result,
                                    std::set<hash_type>& union_nodes) {
            return get_decision_neighbors(graph, this, result, union_nodes);
        }

        bool get_decision_neighbors(GraphType *          graph,
                                    dBG *                shifter,
                                    const std::string&   root,
                                    NeighborBundle&      result,
                                    std::set<hash_type>& union_nodes) {
            shifter->set_cursor(root);
            return get_decision_neighbors(graph, shifter, result, union_nodes);
        }

        bool get_decision_neighbors(GraphType *          graph,
                                    dBG *                shifter,
                                    NeighborBundle&      result,
                                    std::set<hash_type>& union_nodes) {

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
