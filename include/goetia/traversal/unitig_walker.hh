/**
 * (c) Camille Scott, 2019
 * File   : assembly.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 23.07.2019
 */


#ifndef GOETIA_UNITIG_WALKER_HH
#define GOETIA_UNITIG_WALKER_HH

#include "goetia/goetia.hh"
#include "goetia/hashing/kmeriterator.hh"
#include "goetia/hashing/hashextender.hh"
#include "goetia/hashing/canonical.hh"
#include "goetia/storage/storage.hh"

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


namespace goetia {

typedef std::deque<char> Path;
typedef std::vector<std::string> StringVector;


/**
 * @Synopsis  States that result from attempting to take a
 *            step from a given node.
 *
 * STOP_FWD: there are no neighbors in this direction.
 * DECISION_FWD: there is more than one neighbor.
 * STOP_SEEN: there is a single neighbor but it has been seen.
 * DECISION_BKW: There is a single neighbor, but it is a decision in the other direction.
 * STOP_CALLBACK: There is a single neighbor, but it was rejected
 *                by the callback function.
 * BAD_SEED: The node you tried to start at does not exist.
 * GRAPH_ERROR: The graph is structural unsound (basically a panic).
 * STEP: None of the above: we can move.
 *
 */
enum class TraversalState {
    STOP_FWD,
    DECISION_FWD,
    DECISION_BKW,
    STOP_SEEN,
    STOP_MASKED,
    STOP_CALLBACK,
    BAD_SEED,
    GRAPH_ERROR,
    STEP
};


__attribute__((visibility("default")))
inline std::ostream&
operator<<(std::ostream& os, const TraversalState& _this) {
    os << "<TraversalState state=";
    switch(_this) {
        case TraversalState::STOP_FWD:
            os << "STOP_FWD";
            break;
        case TraversalState::DECISION_FWD:
            os << "DECISION_FWD";
            break;
        case TraversalState::DECISION_BKW:
            os << "DECISION_BKW";
            break;
        case TraversalState::STOP_SEEN:
            os << "STOP_SEEN";
            break;
        case TraversalState::STOP_MASKED:
            os << "STOP_MASKED";
            break;
        case TraversalState::STOP_CALLBACK:
            os << "STOP_CALLBACK";
            break;
        case TraversalState::BAD_SEED:
            os << "BAD_SEED";
            break;
        case TraversalState::GRAPH_ERROR:
            os << "GRAPH_ERROR";
            break;
        case TraversalState::STEP:
            os << "STEP";
            break;
        default:
            os << "UNDEFINED";
            break;
    }
    os << ">";
    return os;
}


template<typename HashType>
struct NullWalkFunctor {
    bool operator()(HashType& node) { return true; }
    virtual ~NullWalkFunctor() {}
};

//
// Functor for `walk` that will mask only mid-unitig --
// that is, the walk will be stopped on a node if it
// is in the given set.
//
template<typename HashType>
struct WalkStopper {

    std::set<HashType>& stoppers;

    explicit WalkStopper(std::set<HashType>& stoppers)
        : stoppers(stoppers) {}

    bool operator()(HashType& node) { 
        if (stoppers.count(node)) {
            found = node;
            return false;
        }
        return true;
    }

    std::optional<HashType> get() {
        return found;
    }

    protected:
        std::optional<HashType> found;
};


template<typename HashType>
struct SeenRemover {
    std::set<HashType>& unseen;
    explicit SeenRemover(std::set<HashType>& unseen)
        : unseen(unseen) {}
    bool operator()(HashType& node) {
        unseen.erase(node);
        return true;
    }

};


template <class T>
struct UnitigWalker;


template <template <class, class, class...> class GraphType,
                    class StorageType,
                           class ShifterType,
                                  class... Extras> // Dummy param to handle query in derived classes
class UnitigWalker<GraphType<StorageType, ShifterType, Extras...>> : public extender_selector<ShifterType>::type {

    typedef GraphType<StorageType, ShifterType, Extras...>     Derived;

public:

    typedef Derived                                 graph_type;
    typedef ShifterType                             shifter_type;
    typedef typename extender_selector<shifter_type>::type extender_type;

    typedef typename extender_type::hash_type       hash_type;
    typedef typename hash_type::value_type          value_type;

    template<bool Dir>
        using shift_type = Shift<hash_type, Dir>;
    typedef typename extender_type::shift_left_type  shift_left_type;
    typedef typename extender_type::shift_right_type shift_right_type;

    typedef typename extender_type::kmer_type        kmer_type;

    typedef std::pair<std::vector<kmer_type>,
                      std::vector<kmer_type>>        neighbor_pair_type;

    typedef std::pair<std::vector<shift_type<DIR_LEFT>>,
                      std::vector<shift_type<DIR_RIGHT>>> shift_pair_type;

    typedef TraversalState                    State;
    typedef NullWalkFunctor<hash_type>               null_walk_func_t;

private:

    // dummy template for overload specialization
    // for templated member functions
    template<typename T> struct type { };

public:

    template<bool Dir>
    struct WalkBase {
        kmer_type                                        start;
        std::vector<Shift<hash_type, Dir>> path;
        TraversalState                            end_state;

        const hash_type head() const {
            return this->start;
        }

        const hash_type tail() const {
            if (path.size()) {
                return path.back();
            } else {
                return start;
            }
        }
    };

    template<bool Dir = DIR_RIGHT, typename Dummy = void>
    struct WalkImpl : public WalkBase<Dir> {
        using WalkBase<Dir>::start;
        using WalkBase<Dir>::path;
        using WalkBase<Dir>::end_state;

        const std::string to_string() const {
            std::string str(0, ' ');
            str.reserve(start.kmer.length() + path.size());
            str.assign(start.kmer.begin(), start.kmer.end());
            for (auto& shift : path) {
                str.push_back(shift.symbol);
            }
            return std::move(str);
        }

        const std::string glue(const WalkImpl<DIR_LEFT>& left) const {
            return left.to_string() + this->to_string().substr(start.kmer.size());
        }

        const std::string glue(const WalkImpl<DIR_RIGHT>& right) const {
            return this->to_string() + right.to_string().substr(start.kmer.size());
        }

        const char retreat_symbol() const {
            size_t position = path.size();
            if (position < start.kmer.size()) {
                return start.kmer[position];
            } else {
                return path[position - start.kmer.size()].symbol;
            }
        }
    };

    template<typename Dummy>
    struct WalkImpl<DIR_LEFT, Dummy> : public WalkBase<DIR_LEFT> {
        using WalkBase<DIR_LEFT>::start;
        using WalkBase<DIR_LEFT>::path;
        using WalkBase<DIR_LEFT>::end_state;

        const char retreat_symbol() const {
            size_t position = path.size();
            if (position < start.kmer.size()) {
                return start.kmer[start.kmer.size() - position - 1];
            } else {
                return path[position - start.kmer.size()].symbol;
            }
        }   

        const std::string to_string() const {
            std::string str(0, ' ');
            str.reserve(start.kmer.length() + path.size());

            for (auto it = path.crbegin(); it != path.crend(); ++it) {
                str.push_back(it->symbol);
            }
            std::copy(start.kmer.begin(), start.kmer.end(), std::back_inserter(str));
            return std::move(str);
        }

        const std::string glue(const WalkImpl<DIR_RIGHT>& right) const {
            return this->to_string() + right.to_string().substr(start.kmer.size());
        }
    };

    template<bool Dir>
    struct Walk : public WalkImpl<Dir> {
        typedef WalkImpl<Dir> walk_type;

        using walk_type::start;
        using walk_type::path;
        using walk_type::end_state;
        using walk_type::to_string;
        using walk_type::head;
        using walk_type::tail;
        using walk_type::retreat_symbol;
        using walk_type::glue;
    };

    typedef std::pair<Walk<DIR_LEFT>, Walk<DIR_RIGHT>> walk_pair_type;


    using extender_type::set_cursor;
    using extender_type::get_cursor;
    using extender_type::get;
    using extender_type::hash;
    using extender_type::left_extensions;
    using extender_type::right_extensions;
    using extender_type::shift_left;
    using extender_type::shift_right;
    using extender_type::K;

    std::set<value_type> seen;

    void clear_seen() {
        seen.clear();
    }

    size_t in_degree() {
        auto extensions = this->left_extensions();
        return count_nodes(extensions);
    }

    size_t out_degree() {
        auto extensions = this->right_extensions();
        return count_nodes(extensions);
    }

    size_t degree() {
        return in_degree() + out_degree();
    }

    size_t in_degree(std::set<hash_type>& extras) {

        auto extensions = this->left_extensions();
        return count_nodes(extensions, extras);
    }

    size_t out_degree(std::set<hash_type>& extras) {

        auto extensions = this->right_extensions();
        return count_nodes(extensions, extras);
    }

    size_t degree(std::set<hash_type>& extras) {
        return in_degree(extras) + out_degree(extras);
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
    template<bool Dir>
    size_t count_nodes(const std::vector<shift_type<Dir>>& extensions) {

        uint8_t n_found = 0;
        for (const auto& ext : extensions) {
            if(derived().query(ext)) {
                ++n_found;
            }
        }
        return n_found;
    }


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
    template<bool Dir>
    size_t count_nodes(const std::vector<shift_type<Dir>>& extensions,
                       std::set<hash_type>&                extras) {

        uint8_t n_found = 0;
        for (const auto& ext : extensions) {
            if(derived().query(ext) ||
               extras.count(ext)) {
                ++n_found;
            }
        }
        return n_found;
    }

    /**
     * @Synopsis  Return only the shifts from nodes that exist in the graph.
     *
     * @Param graph
     * @Param nodes
     *
     * @Returns   The valid shifts.
     */
    template<bool Dir>
    auto filter_nodes(const std::vector<shift_type<Dir>>& extensions)
    -> std::vector<shift_type<Dir>> {

        std::vector<shift_type<Dir>> result;
        for (const auto& ext : extensions) {
            if (derived().query(ext)) {
                result.push_back(ext);
            }
        }
        return result;
    }

    /**
     * @Synopsis  Return only the shifts from nodes that exist in the graph.
     *
     * @Param graph
     * @Param nodes
     *
     * @Returns   The valid shifts.
     */
    template<bool Dir>
    auto filter_nodes(const std::vector<shift_type<Dir>>& extensions,
                      std::vector<count_t>& counts)
    -> std::vector<shift_type<Dir>> {

        std::vector<shift_type<Dir>> result;
        for (const auto& ext : extensions) {
            auto count = derived().query(ext);
            if (count != 0) {
                result.push_back(ext);
                counts.push_back(count);
            }
        }
        return result;
    }



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
    template<bool Dir>
    auto filter_nodes(const std::vector<shift_type<Dir>>& extensions,
                      std::set<hash_type>&                extra)
    -> std::vector<shift_type<Dir>> {
        
        std::vector<shift_type<Dir>> result;
        for (const auto& ext : extensions) {
            if (derived().query(ext) ||
                extra.count(ext)) {
                result.push_back(ext);
            }
        }
        return result;
    }


    /**
     * @Synopsis  Filter from a pair of shift vectors (convenience for traversal)
     *
     * @Param graph
     * @Param nodes
     *
     * @Returns   
     */
    shift_pair_type filter_nodes(const shift_pair_type& extensions) {
        return std::make_pair(filter_nodes(extensions.first),
                              filter_nodes(extensions.second));
    }

    shift_pair_type filter_nodes(const shift_pair_type& extensions,
                                 std::set<hash_type>&   extras) {
        return std::make_pair(filter_nodes(extensions.first, extras),
                              filter_nodes(extensions.second, extras));
    }

    /**
     * @Synopsis  In-neighbors for the node given by the current cursor.
     *
     * @Param graph
     *
     * @Returns   Vector of valid shifts in the graph.
     */
    template<class Ret = shift_type<DIR_LEFT>>
    inline auto in_neighbors()
    -> std::vector<Ret> {
        
        return in_neighbors(type<Ret>());
    }

    template<class Ret = shift_type<DIR_LEFT>>
    inline auto in_neighbors(const std::string& seed)
    -> std::vector<Ret> {
        this->set_cursor(seed);
        return in_neighbors(type<Ret>());
    }

private:

    template<class Ret = shift_type<DIR_LEFT>>
    inline auto in_neighbors(type<Ret>)
    -> std::vector<Ret> {
        return filter_nodes(this->left_extensions());
    }

    inline auto in_neighbors(type<kmer_type>)
    -> std::vector<kmer_type> {
        auto nodes = in_neighbors<>();
        return build_left_kmers(nodes, this->get_cursor());
    }

public:

    /**
     * @Synopsis  Out-neighbors for the node given by the current cursor.
     *
     * @Param graph
     *
     * @Returns   Vector of valid shifts in the graph.
     */
    template<class Ret = shift_type<DIR_RIGHT>>
    inline auto out_neighbors()
    -> std::vector<Ret> {
        
        return out_neighbors(type<Ret>());
    }

    template<class Ret = shift_type<DIR_RIGHT>>
    inline auto out_neighbors(const std::string& seed)
    -> std::vector<Ret> {
        this->set_cursor(seed);
        return out_neighbors(type<Ret>());
    }

private:

    template<class Ret = shift_type<DIR_RIGHT>>
    inline auto out_neighbors(type<Ret>)
    -> std::vector<Ret> {
        return filter_nodes(this->right_extensions());
    }

    inline auto out_neighbors(type<kmer_type>)
    -> std::vector<kmer_type> {
        auto nodes = out_neighbors<>();
        return build_right_kmers(nodes, this->get_cursor());
    }

public:


    /**
     * @Synopsis  In and out-neighbors, respectively, for the node
     *            given by the current cursor.
     *
     * @Param graph
     *
     * @Returns   pair of vectors of valid shifts in the graph.
     */
    template<class Ret = shift_pair_type>
    inline Ret neighbors() {
        auto _in = in_neighbors<>();
        auto _out = out_neighbors<>();
        return std::make_pair(_in, _out);
    }

    template<class Ret = shift_pair_type>
    inline Ret neighbors(const std::string& seed) {
        this->set_cursor(seed);
        auto _in = in_neighbors<>();
        auto _out = out_neighbors<>();
        return std::make_pair(_in, _out);
    }

private:

    inline auto neighbors(type<kmer_type>)
    -> neighbor_pair_type {
        auto filtered = neighbors<>();
        auto root = this->get_cursor();
        return std::make_pair(build_left_kmers(filtered.first, root),
                              build_right_kmers(filtered.second, root));
    }



public:

    /**
     * @Synopsis  Given a root k-mers and the shift bases to its left,
     *            build the k-mer strings that could be its left neighbors.
     *
     * @Param nodes The shift_t objects containing the prefix bases and neighbor hashes.
     * @Param root The root k-mer.
     *
     * @Returns   kmer_t objects with the left k-mers.
     */
    std::vector<kmer_type> build_left_kmers(const std::vector<shift_type<DIR_LEFT>>& nodes,
                                            const std::string&                       root) {
        std::vector<kmer_type> kmers;
        auto _prefix = derived().prefix(root);
        for (const auto& neighbor : nodes) {
            kmers.push_back(kmer_type(neighbor,
                                      neighbor.symbol + _prefix));
        }
        return kmers;
    }

    /**
     * @Synopsis  Given a root k-mers and the shift bases to its right,
     *            build the k-mer strings that could be its right neighbors.
     *
     * @Param nodes The shift_t obvjects containing the suffix bases and neighbor hashes.
     * @Param root The root k-mer.
     *
     * @Returns   kmer_t objects with the right k-mers.
     */
    std::vector<kmer_type> build_right_kmers(const std::vector<shift_type<DIR_RIGHT>>& nodes,
                                             const std::string& root) {
        std::vector<kmer_type> kmers;
        auto _suffix = derived().suffix(root);
        for (const auto& neighbor : nodes) {
            kmers.push_back(kmer_type(neighbor,
                                     _suffix + neighbor.symbol));
        }

        return kmers;
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
    template<bool Dir>
    State look_state(std::vector<shift_type<Dir>>& neighbors) {
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
     * @brief Take a single step in DIR_LEFT, if possible. WalkFunctor is a
     *        a callback that takes `hash_type` and returns true if the node
     *        is admitted, and false if not. It is only triggered if STEP_LEFT
     *        is satisfied, meaning it is used to check nodes inside unitigs.
     * 
     * @tparam WalkFunctor 
     * @param f 
     * @return std::pair<State, std::vector<shift_type<DIR_LEFT>>> 
     */
    template<typename WalkFunctor>
    auto step_left(WalkFunctor& f)
    -> std::pair<State, std::vector<shift_type<DIR_LEFT>>> {

        auto neighbors = this->filter_nodes(this->left_extensions());
        auto state = look_state(neighbors);
        if (state == State::STEP) {
            if (!f(neighbors.front().hash)) {
                return {State::STOP_CALLBACK, std::move(neighbors)};
            }

            this->shift_left(neighbors.front().symbol);
            this->seen.insert(neighbors.front().value());
        }
        return {state, std::move(neighbors)};
    }

    /**
     * @brief Overload for step_left with default null functor for WalkFunctor.
     * 
     * @return std::pair<State, std::vector<shift_type<DIR_LEFT>>> 
     */
    auto step_left()
    -> std::pair<State, std::vector<shift_type<DIR_LEFT>>> {
        null_walk_func_t func;
        return step_left(func);
    }

    /**
     * @brief Take a single step in DIR_RIGHT, if possible. WalkFunctor is a
     *        a callback that takes `hash_type` and returns true if the node
     *        is admitted, and false if not. It is only triggered if STEP_LEFT
     *        is satisfied, meaning it is used to check nodes inside unitigs.
     * 
     * @tparam WalkFunctor 
     * @param f 
     * @return std::pair<State, std::vector<shift_type<DIR_RIGHT>>> 
     */
    template<typename WalkFunctor>
    auto step_right(WalkFunctor& f)
    -> std::pair<State, std::vector<shift_type<DIR_RIGHT>>> {

        auto neighbors = this->filter_nodes(this->right_extensions());
        auto state = look_state(neighbors);
        if (state == State::STEP) {
            if (!f(neighbors.front().hash)) {
                return {State::STOP_CALLBACK, std::move(neighbors)};
            }

            this->shift_right(neighbors.front().symbol);
            this->seen.insert(neighbors.front().value());
        }
        return {state, std::move(neighbors)};
    }

    /**
     * @brief Overload for step_right with default null functor for WalkFunctor.
     * 
     * @return std::pair<State, std::vector<shift_type<DIR_RIGHT>>> 
     */
    auto step_right()
    -> std::pair<State, std::vector<shift_type<DIR_RIGHT>>> {
        null_walk_func_t func;
        return step_right(func);
    }

    /**
     * @brief Make as many steps left as possible, starting at seed,
     *        stopping when a STOP state is encountered. If seed
     *        does not exist in the graph, terminates immediately
     *        on BAD_SEED. Stops if WalkFunctor returns false on
     *        a State::STEP.
     * 
     * @tparam WalkFunctor 
     * @param seed 
     * @param f 
     * @return Walk<DIR_LEFT> 
     */
    template<typename WalkFunctor>
    Walk<DIR_LEFT> walk_left(const std::string& seed,
                                      WalkFunctor& f) {

        this->set_cursor(seed);
        auto seed_hash = this->get();
        if (!derived().query(seed_hash)) {
            Walk<DIR_LEFT> walk{kmer_type{seed_hash, this->get_cursor()},
                                         {},
                                         State::BAD_SEED};
            return walk;
        }

        return walk_left(f);
    }

    /**
     * @brief Overload for walk_left with default null walk functor.
     * 
     * @param seed 
     * @return Walk<DIR_LEFT> 
     */
    Walk<DIR_LEFT> walk_left(const std::string& seed) {
        null_walk_func_t func;
        return walk_left(seed, func);
    }

    /**
     * @brief Make as many steps left as possible, starting at cursor,
     *        stopping when a STOP state is encountered. If seed
     *        does not exist in the graph, terminates immediately
     *        on BAD_SEED. Stops if WalkFunctor returns false on
     *        a State::STEP.
     * 
     * @tparam WalkFunctor 
     * @param f 
     * @return Walk<DIR_LEFT> 
     */
    template<typename WalkFunctor>
    Walk<DIR_LEFT> walk_left(WalkFunctor& f) {

        Walk<DIR_LEFT> walk;

        hash_type start_hash = this->get();
        this->seen.clear();
        this->seen.insert(start_hash.value());
        walk.start.hash = start_hash;
        walk.start.kmer = this->get_cursor();

        // take the first step without check for reverse d-nodes
        auto [state, neighbors] = step_left(f);

        if (state != State::STEP) {
            walk.end_state = state;
            return walk;
        } else {
            walk.path.push_back(neighbors.front());
        }

        while (1) {
            if (out_degree() > 1) {
                pdebug("Stop: reverse d-node");
                walk.path.pop_back();
                walk.end_state = State::DECISION_BKW;
                this->seen.erase(this->get().value());
                return std::move(walk);
            }

            const auto [state, neighbors] = step_left(f);

            if (state != State::STEP) {
                walk.end_state = state;
                
                return std::move(walk);
            } else {
                walk.path.push_back(neighbors.front());
            }
        }
    }

    /**
     * @brief Overload for walk_left with default null walk functor.
     * 
     * @return Walk<DIR_LEFT> 
     */
    Walk<DIR_LEFT> walk_left() {
        null_walk_func_t func;
        return walk_left(func);
    }

    /**
     * @brief Make as many steps right as possible, starting at seed,
     *        stopping when a STOP state is encountered. If seed
     *        does not exist in the graph, terminates immediately
     *        on BAD_SEED. Stops if WalkFunctor returns false on
     *        a State::STEP.
     * 
     * @tparam WalkFunctor 
     * @param seed 
     * @param f 
     * @return Walk<DIR_RIGHT> 
     */
    template<typename WalkFunctor>
    Walk<DIR_RIGHT> walk_right(const std::string& seed,
                                        WalkFunctor& f) {

        this->set_cursor(seed);
        auto seed_hash = this->get();
        if (!derived().query(seed_hash)) {
            Walk<DIR_RIGHT> walk{kmer_type{seed_hash, this->get_cursor()},
                                           {},
                                           State::BAD_SEED};
            return walk;
        }
        return walk_right(f);
    }

    /**
     * @brief Overload for walk_right with default null walk functor.
     * 
     * @param seed 
     * @return Walk<DIR_RIGHT> 
     */
    Walk<DIR_RIGHT> walk_right(const std::string& seed) {
        null_walk_func_t func;
        return walk_right(seed, func);
    }

    /**
     * @brief Make as many steps right as possible, starting at cursor,
     *        stopping when a STOP state is encountered. If seed
     *        does not exist in the graph, terminates immediately
     *        on BAD_SEED. Stops if WalkFunctor returns false on
     *        a State::STEP.
     * 
     * @tparam WalkFunctor 
     * @param f 
     * @return Walk<DIR_RIGHT> 
     */
    template<typename WalkFunctor>
    Walk<DIR_RIGHT> walk_right(WalkFunctor& f) {

        Walk<DIR_RIGHT> walk;

        hash_type start_hash = this->get();
        this->seen.clear();
        this->seen.insert(start_hash.value());
        walk.start.hash = start_hash;
        walk.start.kmer = this->get_cursor();

        auto step = step_right(f);

        if (step.first != State::STEP) {
            walk.end_state = step.first;
            return walk;
        } else {
            walk.path.push_back(step.second.front());
        }
        
        while (1) {
            if (in_degree() > 1) {
                walk.path.pop_back();
                walk.end_state = State::DECISION_BKW;
                this->seen.erase(this->get().value());
                return std::move(walk);
            }

            const auto [state, neighbors] = step_right(f);

            if (state != State::STEP) {
                walk.end_state = state;
                return std::move(walk);
            } else {
                walk.path.push_back(neighbors.front());
            }
        }
    }

    /**
     * @brief Overload for walk_right with default null walk functor.
     * 
     * @return Walk<DIR_RIGHT> 
     */
    Walk<DIR_RIGHT> walk_right() {
        null_walk_func_t func;
        return walk_right(func);
    }

    /**
     * @brief Make both left and right walks starting at seed.
     * 
     * @tparam WalkFunctor 
     * @param seed 
     * @param f 
     * @return walk_pair_type 
     */
    template<typename WalkFunctor>
    walk_pair_type walk(const std::string& seed,
                        WalkFunctor& f) {

        this->set_cursor(seed);
        auto seed_hash = this->get();
        if (!derived().query(seed_hash)) {
            kmer_type start{seed_hash, this->get_cursor()};
            return {{start, {}, State::BAD_SEED},
                    {start, {}, State::BAD_SEED}};
        }

        auto left_walk = walk_left(f);
        this->set_cursor(seed);
        auto right_walk = walk_right(f);
        return {left_walk, right_walk};
    }

    /**
     * @brief Default overload of walk with null_walk_func_t.
     * 
     * @param seed 
     * @return walk_pair_type 
     */
    walk_pair_type walk(const std::string& seed) {
        null_walk_func_t func;
        return walk(seed, func);
    }

    std::vector<std::string> extract_unitigs(const std::string& sequence,
                                             std::vector<hash_type>& hashes) {
        

        std::vector<std::string> unitigs;
        null_walk_func_t f;
        std::set<value_type> seq_seen;
        for (size_t i = 0; i < sequence.size() - this->K + 1; ++i) {
            if (seq_seen.count(hashes[i].value())) {
                continue;
            }
            auto [lwalk, rwalk] = walk(sequence.substr(i, this->K), f);
            seq_seen.insert(this->seen.begin(), this->seen.end());
            unitigs.push_back(lwalk.glue(rwalk));
        }

        return unitigs;
    }
                  
    bool is_decision_kmer(const std::string& node,
                          uint8_t&           degree) {

        this->set_cursor(node);
        return is_decision_kmer(degree);
    }

    bool is_decision_kmer(const std::string& node) {

        this->set_cursor(node);
        return this->in_degree() > 1 || this->out_degree() > 1;
    }

    bool is_decision_kmer(uint8_t&    degree) {

        uint8_t ldegree, rdegree;
        ldegree = this->in_degree();
        rdegree = this->out_degree();
        degree = ldegree + rdegree;
        return ldegree > 1 || rdegree > 1;
    }

    void find_decision_kmers(const std::string&           sequence,
                             std::vector<uint32_t>&       decision_positions,
                             std::vector<hash_type>&      decision_hashes,
                             std::vector<neighbor_pair_type>& decision_neighbors) {

        KmerIterator<UnitigWalker> iter(sequence, this);
        size_t pos = 0;
        while(!iter.done()) {
            hash_type h = iter.next();
            neighbor_pair_type neighbors;
            if (get_decision_neighbors(iter.shifter,
                                       neighbors)) {

                decision_neighbors.push_back(neighbors);
                decision_positions.push_back(pos);
                decision_hashes.push_back(h);
            }
        
            ++pos;
       }
    }

    bool get_decision_neighbors(const std::string&   root,
                                neighbor_pair_type&      result,
                                std::set<hash_type>& union_nodes) {

        return get_decision_neighbors(this, root, result, union_nodes);
    }

    bool get_decision_neighbors(neighbor_pair_type&      result,
                                std::set<hash_type>& union_nodes) {
        return get_decision_neighbors(this, result, union_nodes);
    }

    bool get_decision_neighbors(extender_type *      extender,
                                const std::string&   root,
                                neighbor_pair_type&  result,
                                std::set<hash_type>& union_nodes) {
        extender->set_cursor(root);
        return get_decision_neighbors(extender, result, union_nodes);
    }

    bool get_decision_neighbors(extender_type *      extender,
                                neighbor_pair_type&  result,
                                std::set<hash_type>& union_nodes) {
        auto root = extender->get_cursor();
        auto _lfiltered = filter_nodes(extender->left_extensions(), union_nodes);
        auto left_kmers = this->build_left_kmers(_lfiltered, root);
         
        auto _rfiltered = filter_nodes(extender->right_extensions(), union_nodes);
        auto right_kmers = this->build_right_kmers(_rfiltered, root);
        
        if (left_kmers.size() > 1 || right_kmers.size() > 1) {
            result = std::make_pair(left_kmers, right_kmers);
            return true;
        } else {
            return false;
        }
    }

    bool get_decision_neighbors(const std::string& root,
                                neighbor_pair_type&    result) {
        return get_decision_neighbors(this, root, result);
    }

    bool get_decision_neighbors(neighbor_pair_type& result) {
        return get_decision_neighbors(this, result);
    }

    bool get_decision_neighbors(extender_type *     extender,
                                const std::string&  root,
                                neighbor_pair_type& result) {
        extender->set_cursor(root);
        return get_decision_neighbors(extender, result);
    }

    bool get_decision_neighbors(extender_type *     extender,
                                neighbor_pair_type& result) {

        auto root = extender->get_cursor();
        auto _lfiltered = filter_nodes(extender->left_extensions());
        auto left_kmers = this->build_left_kmers(_lfiltered, root);
         
        auto _rfiltered = filter_nodes(extender->right_extensions());
        auto right_kmers = this->build_right_kmers(_rfiltered, root);

        if (left_kmers.size() > 1 || right_kmers.size() > 1) {
            result = std::make_pair(left_kmers, right_kmers);
            return true;
        } else {
            return false;
        }
    }

protected:

    template<typename... ExtraArgs>
    explicit UnitigWalker(uint16_t K, ExtraArgs&&... args)
        : extender_type(K, std::forward<ExtraArgs>(args)...)
    {
    }

    template<typename... ExtraArgs>
    explicit UnitigWalker(const std::string& start,
                       uint16_t K,
                       ExtraArgs&&... args)
        : extender_type(start, K, std::forward<ExtraArgs>(args)...)
    {
    }

    explicit UnitigWalker(const extender_type& extender)
        : extender_type(extender)
    {
    }

    explicit UnitigWalker(const shifter_type& shifter)
        : extender_type(shifter)
    {
    }

    friend Derived;

    Derived& derived() {
        return *static_cast<Derived*>(this);
    }

    const Derived& derived() const {
        return *static_cast<const Derived*>(this);
    }

};


}

#undef pdebug
#endif
