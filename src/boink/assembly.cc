#include "boink/assembly.hh"

#include "boink/dbg.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/sparseppstorage.hh"

# ifdef DEBUG_ASMLY
#   define pdebug(x) do { std::cerr << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif

namespace boink {

template <class GraphType>
size_t
Traverse<GraphType>::dBG::count_nodes(GraphType *                 graph,
                           const std::vector<shift_t>& nodes) {
    uint8_t n_found = 0;
    for (auto node: nodes) {
        if(graph->query(node.hash)) {
            ++n_found;
        }
    }
    return n_found;
}


template <class GraphType>
size_t
Traverse<GraphType>::dBG::count_nodes(GraphType *                 graph,
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


template <class GraphType>
size_t
Traverse<GraphType>::dBG::reduce_nodes(GraphType *                 graph,
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

template <class GraphType>
size_t
Traverse<GraphType>::dBG::reduce_nodes(GraphType *                 graph,
                                       const std::vector<shift_t>& nodes,
                                       shift_t&                    result,
                                       std::set<hash_t>&           extra) {
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


template <class GraphType>
std::vector<shift_t> 
Traverse<GraphType>::dBG::filter_nodes(GraphType * graph,
                                       const std::vector<shift_t>& nodes) {
    std::vector<shift_t> result;
    for (auto node : nodes) {
        if (graph->query(node.hash)) {
            result.push_back(node);
        }
    }
    return result;
}


template <class GraphType>
std::vector<shift_t> 
Traverse<GraphType>::dBG::filter_nodes(GraphType *                 graph,
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


template <class GraphType>
typename Traverse<GraphType>::EndState
Traverse<GraphType>::dBG::traverse_left(GraphType *        graph,
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


template <class GraphType>
typename Traverse<GraphType>::EndState
Traverse<GraphType>::dBG::traverse_left(GraphType * graph,
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


template <class GraphType>
typename Traverse<GraphType>::EndState
Traverse<GraphType>::dBG::traverse_right(GraphType *        graph,
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


template <class GraphType>
typename Traverse<GraphType>::EndState
Traverse<GraphType>::dBG::traverse_right(GraphType *       graph,
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


template <class GraphType>
std::pair<typename Traverse<GraphType>::EndState, typename Traverse<GraphType>::EndState>
Traverse<GraphType>::dBG::traverse(GraphType *        graph,
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

} // namespace boink

template struct boink::Traverse<boink::dBG<boink::storage::BitStorage,
                                           boink::hashing::RollingHashShifter>>;
template struct boink::Traverse<boink::dBG<boink::storage::ByteStorage,
                                           boink::hashing::RollingHashShifter>>;
template struct boink::Traverse<boink::dBG<boink::storage::NibbleStorage,
                                           boink::hashing::RollingHashShifter>>;
template struct boink::Traverse<boink::dBG<boink::storage::QFStorage,
                                           boink::hashing::RollingHashShifter>>;
template struct boink::Traverse<boink::dBG<boink::storage::SparseppSetStorage,
                                           boink::hashing::RollingHashShifter>>;

#undef pdebug
