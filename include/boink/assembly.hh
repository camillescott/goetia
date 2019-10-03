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
        GRAPH_ERROR
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

        template<typename... Args>
        explicit dBG(Args&&... args)
            : BaseShifter(std::forward<Args>(args)...)
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

        size_t degree_left(GraphType *          graph,
                           std::set<hash_type>& extras) {

            auto neighbors = this->gather_left();
            return count_nodes(graph, neighbors, extras);
        }

        size_t degree_right(GraphType *          graph,
                            std::set<hash_type>& extras) {

            auto neighbors = this->gather_right();
            return count_nodes(graph, neighbors, extras);
        }

        size_t degree(GraphType *          graph,
                      std::set<hash_type>& extras) {
            return degree_left(graph, extras) + degree_right(graph, extras);
        }

        size_t try_traverse_left(GraphType * graph,
                                 shift_type& result) {

            std::vector<shift_type> neighbors = this->gather_left();
            auto n_left = reduce_nodes(graph, neighbors, result);

            if (n_left == 1) {
                this->shift_left(result.symbol);
            }
            return n_left;
        }

        size_t try_traverse_right(GraphType * graph,
                                  shift_type& result) {
            std::vector<shift_type> neighbors = this->gather_right();
            auto n_right = reduce_nodes(graph, neighbors, result);

            if (n_right == 1) {
                this->shift_right(result.symbol);
            }
            return n_right;
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

        EndState traverse_left(GraphType *          graph,
                               const std::string&   seed,
                               Path&                path,
                               std::set<hash_type>& mask);

        EndState traverse_left(GraphType *          graph,
                               Path&                path,
                               std::set<hash_type>& mask);

        EndState traverse_right(GraphType *          graph,
                                const std::string&   seed,
                                Path&                path,
                                std::set<hash_type>& mask);

        EndState traverse_right(GraphType *       graph,
                                Path&             path,
                                std::set<hash_type>& mask);

        std::pair<EndState, EndState> traverse(GraphType *         graph,
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
