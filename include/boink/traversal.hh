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

template <class T>
struct dBGWalker;


template <template <class, class> class GraphType,
                                  class StorageType,
                                  class ShifterType>
class dBGWalker<GraphType<StorageType, ShifterType>> : public hashing::HashExtender<ShifterType> {

    typedef GraphType<StorageType, ShifterType>     Derived;

public:

    typedef Derived                                 graph_type;
    typedef ShifterType                             shifter_type;
    typedef hashing::HashExtender<shifter_type>     extender_type;

    typedef typename extender_type::hash_type       hash_type;
    typedef typename hash_type::value_type          value_type;

    template<bool Dir>
        using shift_type = hashing::ShiftModel<hash_type, Dir>;
    typedef typename extender_type::shift_left_type  shift_left_type;
    typedef typename extender_type::shift_right_type shift_right_type;

    typedef typename extender_type::kmer_type        kmer_type;

    typedef std::pair<std::vector<kmer_type>,
                      std::vector<kmer_type>>        neighbor_pair_type;

    typedef std::pair<std::vector<shift_type<hashing::DIR_LEFT>>,
                      std::vector<shift_type<hashing::DIR_RIGHT>>> shift_pair_type;

    typedef TraversalState::State       State;

private:

    // dummy template for overload specialization
    // for templated member functions
    template<typename T> struct type { };

public:

    template<bool Dir>
    struct WalkBase {
        kmer_type                                        start;
        std::vector<hashing::ShiftModel<hash_type, Dir>> path;
        TraversalState::State                            end_state;

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

    template<bool Dir = hashing::DIR_RIGHT, typename Dummy = void>
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

        const std::string glue(const WalkImpl<hashing::DIR_LEFT>& left) const {
            return left.to_string() + this->to_string().substr(start.kmer.size());
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
            return std::move(str);
        }

        const std::string glue(const WalkImpl<hashing::DIR_RIGHT>& right) const {
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
        using walk_type::glue;
    };

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
     * @Synopsis  Count nodes that short-circuits when more
     *            than one is found (for traversal).
     *
     * @Param graph
     * @Param nodes
     * @Param result The shift if n_bound is one.
     *
     * @Returns   
     *
    template<bool Dir>
    shift_type<Dir> reduce_nodes(GraphType *                         graph,
                                 const std::vector<shift_type<Dir>>& extensions) {

        uint8_t n_found = 0;
        shift_type<Dir> result;
        for (const auto& ext : extensions) {
            if(graph->query(ext)) {
                ++n_found;
                if (n_found > 1) {
                    return n_found;
                }
                result = ext;
            }
        }
        return result;
    }

    template<bool Dir>
    shift_type<Dir> reduce_nodes(GraphType *                         graph,
                                 const std::vector<shift_type<Dir>>& extensions,
                                 std::set<hash_type>&                extra) {
        uint8_t n_found = 0;
        shift_type<Dir> result;
        for (const auto& ext : extensions ) {
            if(graph->query(ext) ||
               extra.count(ext)) {
                ++n_found;
                if (n_found > 1) {
                    return n_found;
                }
                result = ext;
            }
        }
        return result;
    }
    */


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
    template<class Ret = shift_type<hashing::DIR_LEFT>>
    inline auto in_neighbors()
    -> std::vector<Ret> {
        
        return in_neighbors(type<Ret>());
    }

    template<class Ret = shift_type<hashing::DIR_LEFT>>
    inline auto in_neighbors(const std::string& seed)
    -> std::vector<Ret> {
        this->set_cursor(seed);
        return in_neighbors(type<Ret>());
    }

private:

    template<class Ret = shift_type<hashing::DIR_LEFT>>
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
    template<class Ret = shift_type<hashing::DIR_RIGHT>>
    inline auto out_neighbors()
    -> std::vector<Ret> {
        
        return out_neighbors(type<Ret>());
    }

    template<class Ret = shift_type<hashing::DIR_RIGHT>>
    inline auto out_neighbors(const std::string& seed)
    -> std::vector<Ret> {
        this->set_cursor(seed);
        return out_neighbors(type<Ret>());
    }

private:

    template<class Ret = shift_type<hashing::DIR_RIGHT>>
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
    std::vector<kmer_type> build_left_kmers(const std::vector<shift_type<hashing::DIR_LEFT>>& nodes,
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
    std::vector<kmer_type> build_right_kmers(const std::vector<shift_type<hashing::DIR_RIGHT>>& nodes,
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
     * @Synopsis  Gather and filter left neighbors from the current cursor,
     *            then step left if the look_state was STEP.
     *
     * @Param graph
     *
     * @Returns   Pair of the traversal state and any found neighbors.
     */
    auto step_left()
    -> std::pair<State, std::vector<shift_type<hashing::DIR_LEFT>>> {

        auto neighbors = this->filter_nodes(this->left_extensions());
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
    auto step_left(std::set<hash_type>& mask)
    -> std::pair<State, std::vector<shift_type<hashing::DIR_LEFT>>> {

        auto neighbors = this->filter_nodes(this->left_extensions());
        auto state = look_state(neighbors);
        if (state == State::STEP) {
            if (mask.count(neighbors.front())) {
                state = State::STOP_MASKED;
            } else {
                this->shift_left(neighbors.front().symbol);
                this->seen.insert(neighbors.front().value());
            }
        }
        return {state, std::move(neighbors)};
    }

    auto step_right()
    -> std::pair<State, std::vector<shift_type<hashing::DIR_RIGHT>>> {

        auto neighbors = this->filter_nodes(this->right_extensions());
        auto state = look_state(neighbors);
        if (state == State::STEP) {
            this->shift_right(neighbors.front().symbol);
            this->seen.insert(neighbors.front().value());
        }
        return {state, std::move(neighbors)};
    }

    auto step_right(std::set<hash_type>& mask)
    -> std::pair<State, std::vector<shift_type<hashing::DIR_RIGHT>>> {

        auto neighbors = this->filter_nodes(this->right_extensions());
        auto state = look_state(neighbors);
        if (state == State::STEP) {
            if (mask.count(neighbors.front())) {
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
    Walk<hashing::DIR_LEFT> walk_left(const std::string&   seed,
                                      std::set<hash_type>& mask) {

        this->set_cursor(seed);
        auto seed_hash = this->get();
        if (!derived().query(seed_hash)) {
            Walk<hashing::DIR_LEFT> walk{kmer_type{seed_hash, this->get_cursor()},
                                         {},
                                         State::BAD_SEED};
            return walk;
        } 
        return walk_left(mask);
    }


    /**
     * @Synopsis  Walk left, from the current cursor position.
     *
     * @Param graph
     * @Param path
     * @Param mask
     *
     * @Returns   
     */
    Walk<hashing::DIR_LEFT> walk_left(std::set<hash_type>& mask) {

        Walk<hashing::DIR_LEFT> walk;

        hash_type start_hash = this->get();
        this->seen.clear();
        this->seen.insert(start_hash.value());
        walk.start.hash = start_hash;
        walk.start.kmer = this->get_cursor();

        // take the first step without check for reverse d-nodes
        auto step = step_left(mask);

        if (step.first != State::STEP) {
            walk.end_state = step.first;
            return walk;
        } else {
            walk.path.push_back(step.second.front());
        }

        while (1) {
            if (out_degree() > 1) {
                pdebug("Stop: reverse d-node");
                walk.path.pop_back();
                walk.end_state = State::DECISION_RC;
                this->seen.erase(this->get().value());
                return std::move(walk);
            }

            step = step_left(mask);

            if (step.first != State::STEP) {
                walk.end_state = step.first;
                
                return std::move(walk);
            } else {
                walk.path.push_back(step.second.front());
            }
        }
    }

    Walk<hashing::DIR_RIGHT> walk_right(const std::string&   seed,
                                        std::set<hash_type>& mask) {

        this->set_cursor(seed);
        auto seed_hash = this->get();
        if (!derived().query(seed_hash)) {
            Walk<hashing::DIR_RIGHT> walk{kmer_type{seed_hash, this->get_cursor()},
                                           {},
                                           State::BAD_SEED};
            return walk;
        }
        return walk_right(mask);
    }

    Walk<hashing::DIR_RIGHT> walk_right(std::set<hash_type>& mask) {

        Walk<hashing::DIR_RIGHT> walk;

        hash_type start_hash = this->get();
        this->seen.clear();
        this->seen.insert(start_hash.value());
        walk.start.hash = start_hash;
        walk.start.kmer = this->get_cursor();

        auto step = step_right(mask);

        if (step.first != State::STEP) {
            walk.end_state = step.first;
            return walk;
        } else {
            walk.path.push_back(step.second.front());
        }
        
        while (1) {
            if (in_degree() > 1) {
                walk.path.pop_back();
                walk.end_state = State::DECISION_RC;
                this->seen.erase(this->get().value());
                return std::move(walk);
            }

            step = step_right(mask);

            if (step.first != State::STEP) {
                walk.end_state = step.first;
                return std::move(walk);
            } else {
                walk.path.push_back(step.second.front());
            }
        }
    }

    walk_pair_type walk(const std::string&  seed,
                        std::set<hash_type> mask) {

        this->set_cursor(seed);
        auto seed_hash = this->get();
        if (!derived().query(seed_hash)) {
            kmer_type start{seed_hash, this->get_cursor()};
            return {{start, {}, State::BAD_SEED},
                    {start, {}, State::BAD_SEED}};
        }

        auto left_walk = walk_left(mask);
        this->set_cursor(seed);
        auto right_walk = walk_right(mask);
        return {left_walk, right_walk};
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

        hashing::KmerIterator<dBGWalker> iter(sequence, this);
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

    explicit dBGWalker(const shifter_type& shifter)
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
