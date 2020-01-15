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
#include "boink/hashing/hashextender.hh"
#include "boink/hashing/canonical.hh"

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
using hashing::Direction_t;


/**
 * @Synopsis  States that result from attempting to take a
 *            step from a given node.
 *
 * STOP_FWD: there are no neighbors in this direction.
 * DECISION_FWD: there is more than one neighbor.
 * STOP_SEEN: there is a single neighbor but it has been seen.
 * DECISION_RC: There is a single neighbor, but it is a decision in the other direction.
 * STOP_MASKED: There is a single neighbor, but it is masked.
 * BAD_SEED: The node you tried to start at does not exist.
 * GRAPH_ERROR: The graph is structural unsound (basically a panic).
 * STEP: None of the above: we can move.
 *
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
class dBGWalker : public hashing::HashExtender<typename GraphType::shifter_type> {

public:

    typedef GraphType                               graph_type;
    typedef typename graph_type::shifter_type       shifter_type;
    typedef hashing::HashExtender<shifter_type>     extender_type;

    typedef typename extender_type::alphabet        alphabet;
    typedef typename extender_type::hash_type       hash_type;
    typedef typename hash_type::value_type          value_type;

    template<Direction_t D>
        using shift_type = hashing::ShiftModel<hash_type, D>;
    typedef typename extender_type::shift_left_type  shift_left_type;
    typedef typename extender_type::shift_right_type shift_right_type;

    typedef typename extender_type::kmer_type        kmer_type;

    typedef std::pair<std::vector<kmer_type>,
                      std::vector<kmer_type>>                      neighbor_pair_type;
    typedef std::pair<std::vector<shift_type<hashing::DIR_LEFT>>,
                      std::vector<shift_type<hashing::DIR_RIGHT>>> shift_pair_type;

    typedef TraversalState::State       State;

protected:

    template<Direction_t D>
    struct WalkBase {
        kmer_type                  start;
        std::vector<shift_type<D>> path;
        State                      end_state;

        const hash_type head() const {
            return this->start.value();
        }

        const hash_type tail() const {
            return path.back().value();
        }

    };

    template<Direction_t D = hashing::DIR_RIGHT, typename Dummy = void>
    struct WalkImpl : public WalkBase<D> {
        using WalkBase<D>::start;
        using WalkBase<D>::path;
        using WalkBase<D>::end_state;

        const std::string to_string() const {
            std::string str(0, ' ');
            str.reserve(start.kmer.length() + path.size());
            str.assign(start.kmer.begin(), start.kmer.end());
            for (auto& shift : path) {
                str.push_back(shift.symbol);
            }
            return std::move(str);
        }
    };

    template<typename Dummy>
    struct WalkImpl<hashing::DIR_LEFT, Dummy> : public WalkBase<hashing::DIR_LEFT> {
        using WalkBase<hashing::DIR_LEFT>::start;
        using WalkBase<hashing::DIR_LEFT>::path;
        using WalkBase<hashing::DIR_LEFT>::end_state;

        const std::string to_string() const {
            std::string str(0, ' ');
            str.reserve(start.kmer.length() + path.size());

            for (auto it = path.crbegin(); it != path.crend(); ++it) {
                str.push_back(it->symbol);
            }
            std::copy(start.kmer.begin(), start.kmer.end(), std::back_inserter(str));
        }
    };


    std::set<hash_type> seen;

    using extender_type::_K;

public:

    template<Direction_t D>
        using Walk = WalkImpl<D>;
    typedef std::pair<Walk<hashing::DIR_LEFT>, Walk<hashing::DIR_RIGHT>> walk_pair_type;

    using extender_type::set_cursor;
    using extender_type::get_cursor;
    using extender_type::get;
    using extender_type::hash;
    using extender_type::left_extensions;
    using extender_type::right_extensions;
    using extender_type::shift_left;
    using extender_type::shift_right;
    using extender_type::K;

    template<typename... ExtraArgs>
    explicit dBGWalker(uint16_t K, ExtraArgs&&... args)
        : extender_type(K, std::forward<ExtraArgs>(args)...)
    {
    }

    template<typename... ExtraArgs>
    explicit dBGWalker(const std::string& start,
                       uint16_t K,
                       ExtraArgs&&... args)
        : extender_type(start, K, std::forward<ExtraArgs>(args)...)
    {
    }

    explicit dBGWalker(const extender_type& extender)
        : extender_type(extender)
    {
    }

    template<typename... Args>
    static std::shared_ptr<dBGWalker> build(Args&&... args) {
        return std::make_shared<dBGWalker>(std::forward<Args>(args)...);
    }

    static std::shared_ptr<dBGWalker> build(const extender_type& extender) {
        return std::make_shared<dBGWalker>(extender);
    }

    void clear_seen() {
        seen.clear();
    }

    size_t in_degree(GraphType * graph) {
        auto extensions = this->left_extensions();
        return count_nodes(graph, extensions);
    }

    size_t out_degree(GraphType * graph) {
        auto extensions = this->right_extensions();
        return count_nodes(graph, extensions);
    }

    size_t degree(GraphType * graph) {
        return in_degree(graph) + out_degree(graph);
    }

    size_t in_degree(GraphType *          graph,
                     std::set<hash_type>& extras) {

        auto extensions = this->left_extensions();
        return count_nodes(graph, extensions, extras);
    }

    size_t out_degree(GraphType *          graph,
                      std::set<hash_type>& extras) {

        auto extensions = this->right_extensions();
        return count_nodes(graph, extensions, extras);
    }

    size_t degree(GraphType *          graph,
                  std::set<hash_type>& extras) {
        return in_degree(graph, extras) + out_degree(graph, extras);
    }

    /**
     * @Synopsis  Count how many nodes from the vector
     *            exist in the graph.
     *
     * @Param graph
     * @Param nodes
     *
     * @Returns   Number of nodes that exist in graph.
     */
    template<Direction_t D>
    size_t count_nodes(GraphType *                       graph,
                       const std::vector<shift_type<D>>& extensions);


    /**
     * @Synopsis  Count how many nodes are in the in the induced
     *            graph resulting from the input graph and the 
     *            nodes and implicit edges from extras.
     *
     * @Param graph
     * @Param nodes
     * @Param extras
     *
     * @Returns   Number of nodes.
     */
    template<Direction_t D>
    size_t count_nodes(GraphType *                       graph,
                       const std::vector<shift_type<D>>& extensions,
                       std::set<hash_type>&              extras);

    /**
     * @Synopsis  Count nodes that short-circuits when more
     *            than one is found (for traversal).
     *
     * @Param graph
     * @Param nodes
     * @Param result The shift if n_bound is one.
     *
     * @Returns   
     */
    template<Direction_t D>
    shift_type<D> reduce_nodes(GraphType *                       graph,
                               const std::vector<shift_type<D>>& extensions);

    template<Direction_t D>
    shift_type<D> reduce_nodes(GraphType *                       graph,
                               const std::vector<shift_type<D>>& extensions,
                               std::set<hash_type>&              extra);


    /**
     * @Synopsis  Return only the shifts from nodes that exist in the graph.
     *
     * @Param graph
     * @Param nodes
     *
     * @Returns   The valid shifts.
     */
    template<Direction_t D>
    auto filter_nodes(GraphType *                       graph,
                      const std::vector<shift_type<D>>& extensions)
    -> std::vector<shift_type<D>>;


    /**
     * @Synopsis  Return only the shifts from nodes that exist in the induced
     *            graph of extra on to graph.
     *
     * @Param graph
     * @Param nodes
     * @Param extra
     *
     * @Returns   
     */
    template<Direction_t D>
    auto filter_nodes(GraphType *                       graph,
                       const std::vector<shift_type<D>>& extensions,
                       std::set<hash_type>&              extra)
    -> std::vector<shift_type<D>>;


    /**
     * @Synopsis  Filter from a pair of shift vectors (convenience for traversal)
     *
     * @Param graph
     * @Param nodes
     *
     * @Returns   
     */
    shift_pair_type filter_nodes(GraphType *            graph,
                                 const shift_pair_type& nodes);

    shift_pair_type filter_nodes(GraphType *            graph,
                                 const shift_pair_type& nodes,
                                 std::set<hash_type>&   extras);

    /**
     * @Synopsis  In-neighbors for the node given by the current cursor.
     *
     * @Param graph
     *
     * @Returns   Vector of valid shifts in the graph.
     */
    auto in_neighbors(GraphType * graph)
    -> std::vector<shift_type<hashing::DIR_LEFT>> {
        return filter_nodes(graph, this->left_extensions());
    }

    /**
     * @Synopsis  Out-neighbors for the node given by the current cursor.
     *
     * @Param graph
     *
     * @Returns   Vector of valid shifts in the graph.
     */
    auto out_neighbors(GraphType * graph)
    -> std::vector<shift_type<hashing::DIR_RIGHT>> {
        return filter_nodes(graph, this->right_extensions());
    }


    /**
     * @Synopsis  In and out-neighbors, respectively, for the node
     *            given by the current cursor.
     *
     * @Param graph
     *
     * @Returns   pair of vectors of valid shifts in the graph.
     */
    shift_pair_type neighbors(GraphType * graph) {
        auto _in = in_neighbors(graph);
        auto _out = out_neighbors(graph);
        return std::make_pair(_in, _out);
    }

    std::vector<kmer_type> find_left_kmers(GraphType * graph) {
        auto root = this->get_cursor();
        auto filtered = filter_nodes(graph, this->left_extensions());
        return graph->build_left_kmers(filtered, root);
    }

    std::vector<kmer_type> find_right_kmers(GraphType * graph) {
        auto root = this->get_cursor();
        auto filtered = filter_nodes(graph, this->right_extensions());
        return graph->build_right_kmers(filtered, root);
    }

    std::vector<kmer_type> find_left_kmers(GraphType *       graph,
                                           std::set<hash_type>& extras) {
        auto root = this->get_cursor();
        auto filtered = filter_nodes(graph, this->left_extensions(), extras);
        return graph->build_left_kmers(filtered, root);
    }

    std::vector<kmer_type> find_right_kmers(GraphType *       graph,
                                            std::set<hash_type>& extras) {
        auto root = this->get_cursor();
        auto filtered = filter_nodes(graph, this->right_extensions(), extras);
        return graph->build_right_kmers(filtered, root);
    }

    /**
     * @Synopsis  Given a vector of neighbors prior to stepping, deduce the traversal
     *            state from them. 
     *
     * @Param neighbors
     *
     * @Returns   DECISION_FWD if there is more than 1 neighbor;
     *            STOP_FWD if there are no neighbors;
     *            STOP_SEEN if the single node is in the seen set;
     *            STEP if there is a single node not in the seen set.
     */
    template<Direction_t D>
    State look_state(std::vector<shift_type<D>>& neighbors) {
        if (neighbors.size() > 1) {
            pdebug("Stop: forward d-node");
            return State::DECISION_FWD;
        } else if (neighbors.size() == 0) {
            pdebug("Stop: no neighbors.");
            return State::STOP_FWD;
        } else if (this->seen.count(neighbors.front().value())) {
            pdebug("Stop: hit seen k-mer.");
            return State::STOP_SEEN;
        } else {
            return State::STEP;
        }
    }

    /**
     * @Synopsis  Gather and filter left neighbors from the current cursor,
     *            then step left if the look_state was STEP.
     *
     * @Param graph
     *
     * @Returns   Pair of the traversal state and any found neighbors.
     */
    auto step_left(GraphType * graph)
    -> std::pair<State, std::vector<shift_type<hashing::DIR_LEFT>>> {
        auto neighbors = this->filter_nodes(graph, this->left_extensions());
        auto state = look_state(neighbors);
        if (state == State::STEP) {
            this->shift_left(neighbors.front().symbol);
            this->seen.insert(neighbors.front().value());
        }
        return {state, std::move(neighbors)};
    }

    /**
     * @Synopsis  Step left, while masking the nodes in mask.
     *
     * @Param graph
     * @Param mask
     *
     * @Returns   
     */
    auto step_left(GraphType * graph,
                   std::set<hash_type>& mask)
    -> std::pair<State, std::vector<shift_type<hashing::DIR_LEFT>>> {

        auto neighbors = this->filter_nodes(graph, this->left_extensions());
        auto state = look_state(neighbors);
        if (state == State::STEP) {
            if (mask.count(neighbors.front().value())) {
                state = State::STOP_MASKED;
            } else {
                this->shift_left(neighbors.front().symbol);
                this->seen.insert(neighbors.front().value());
            }
        }
        return {state, std::move(neighbors)};
    }

    auto step_right(GraphType * graph)
    -> std::pair<State, std::vector<shift_type<hashing::DIR_RIGHT>>> {

        auto neighbors = this->filter_nodes(graph, this->right_extensions());
        auto state = look_state(neighbors);
        if (state == State::STEP) {
            this->shift_right(neighbors.front().symbol);
            this->seen.insert(neighbors.front().value());
        }
        return {state, std::move(neighbors)};
    }

    auto step_right(GraphType *          graph,
                    std::set<hash_type>& mask)
    -> std::pair<State, std::vector<shift_type<hashing::DIR_RIGHT>>> {

        auto neighbors = this->filter_nodes(graph, this->right_extensions());
        auto state = look_state(neighbors);
        if (state == State::STEP) {
            if (mask.count(neighbors.front().value())) {
                state = State::STOP_MASKED;
            } else {
                this->shift_right(neighbors.front().symbol);
                this->seen.insert(neighbors.front().value());
            }
        }
        return {state, std::move(neighbors)};
    }

    /**
     * @Synopsis  Make as many steps left as possible, starting at seed,
     *            stopping when a STOP state is encountered. If seed
     *            does not exist in the graph, terminates immediately
     *            on BAD_SEED.
     *
     * @Param graph
     * @Param seed
     * @Param path The string spelled out by the walk.
     * @Param mask Nodes the mask out of the graph.
     *
     * @Returns   
     */
    Walk<hashing::DIR_LEFT> walk_left(GraphType *          graph,
                                      const std::string&   seed,
                                      std::set<hash_type>& mask);


    /**
     * @Synopsis  Walk left, from the current cursor position.
     *
     * @Param graph
     * @Param path
     * @Param mask
     *
     * @Returns   
     */
    Walk<hashing::DIR_LEFT> walk_left(GraphType *          graph,
                                      std::set<hash_type>& mask);

    Walk<hashing::DIR_RIGHT> walk_right(GraphType *          graph,
                                        const std::string&   seed,
                                        std::set<hash_type>& mask);

    Walk<hashing::DIR_RIGHT> walk_right(GraphType *          graph,
                                        std::set<hash_type>& mask);

    walk_pair_type walk(GraphType *         graph,
                        const std::string&  seed,
                        std::set<hash_type> mask);
                  
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
                             std::vector<neighbor_pair_type>& decision_neighbors) {

        hashing::KmerIterator<dBGWalker> iter(sequence, this);
        size_t pos = 0;
        while(!iter.done()) {
            hash_type h = iter.next();
            neighbor_pair_type neighbors;
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
                                neighbor_pair_type&      result,
                                std::set<hash_type>& union_nodes) {

        return get_decision_neighbors(graph, this, root, result, union_nodes);
    }

    bool get_decision_neighbors(GraphType *          graph,
                                neighbor_pair_type&      result,
                                std::set<hash_type>& union_nodes) {
        return get_decision_neighbors(graph, this, result, union_nodes);
    }

    bool get_decision_neighbors(GraphType *          graph,
                                dBGWalker *                extender,
                                const std::string&   root,
                                neighbor_pair_type&      result,
                                std::set<hash_type>& union_nodes) {
        extender->set_cursor(root);
        return get_decision_neighbors(graph, extender, result, union_nodes);
    }

    bool get_decision_neighbors(GraphType *          graph,
                                dBGWalker *                extender,
                                neighbor_pair_type&      result,
                                std::set<hash_type>& union_nodes) {

        auto left_kmers = extender->find_left_kmers(graph, union_nodes);
        auto right_kmers = extender->find_right_kmers(graph, union_nodes);
        
        if (left_kmers.size() > 1 || right_kmers.size() > 1) {
            result = std::make_pair(left_kmers, right_kmers);
            return true;
        } else {
            return false;
        }
    }

    bool get_decision_neighbors(GraphType *        graph,
                                const std::string& root,
                                neighbor_pair_type&    result) {
        return get_decision_neighbors(graph, this, root, result);
    }

    bool get_decision_neighbors(GraphType *     graph,
                                neighbor_pair_type& result) {
        return get_decision_neighbors(graph, this, result);
    }

    bool get_decision_neighbors(GraphType *        graph,
                                dBGWalker *              extender,
                                const std::string& root,
                                neighbor_pair_type&    result) {
        extender->set_cursor(root);
        return get_decision_neighbors(graph, extender, result);
    }

    bool get_decision_neighbors(GraphType *     graph,
                                dBGWalker *           extender,
                                neighbor_pair_type& result) {

        auto left_kmers = extender->find_left_kmers(graph);
        auto right_kmers = extender->find_right_kmers(graph);
        
        if (left_kmers.size() > 1 || right_kmers.size() > 1) {
            result = std::make_pair(left_kmers, right_kmers);
            return true;
        } else {
            return false;
        }
    }

};
}

#undef pdebug
#endif
