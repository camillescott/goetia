/**
 * (c) Camille Scott, 2019
 * File   : assembly.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */
#include "boink/traversal.hh"

#include "boink/dbg.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/ukhs.hh"
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
                           const std::vector<shift_type>& nodes) {
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
                           const std::vector<shift_type>& nodes,
                           std::set<hash_type>&           extras) {
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
                                       const std::vector<shift_type>& nodes,
                                       shift_type&                    result) {
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
                                       const std::vector<shift_type>& nodes,
                                       shift_type&                    result,
                                       std::set<hash_type>&           extra) {
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
std::vector<typename GraphType::shifter_type::shift_type> 
Traverse<GraphType>::dBG::filter_nodes(GraphType * graph,
                                       const std::vector<shift_type>& nodes) {
    std::vector<shift_type> result;
    for (auto node : nodes) {
        if (graph->query(node.hash)) {
            result.push_back(node);
        }
    }
    return result;
}


template <class GraphType>
std::vector<typename GraphType::shifter_type::shift_type> 
Traverse<GraphType>::dBG::filter_nodes(GraphType *                 graph,
                                       const std::vector<shift_type>& nodes,
                                       std::set<hash_type>&           extra) {
    std::vector<shift_type> result;
    for (auto node : nodes) {
        if (graph->query(node.hash) ||
            extra.count(node.hash)) {
            result.push_back(node);
        }
    }
    return result;
}


template <class GraphType>
std::pair<std::vector<typename GraphType::shifter_type::shift_type>,
          std::vector<typename GraphType::shifter_type::shift_type>>
Traverse<GraphType>::dBG::filter_nodes(GraphType * graph,
                                       const std::pair<std::vector<shift_type>,
                                                       std::vector<shift_type>>& nodes) {
    return std::make_pair(filter_nodes(graph, nodes.first),
                          filter_nodes(graph, nodes.second));
}


template <class GraphType>
std::pair<std::vector<typename GraphType::shifter_type::shift_type>,
          std::vector<typename GraphType::shifter_type::shift_type>>
Traverse<GraphType>::dBG::filter_nodes(GraphType *                               graph,
                                       const std::pair<std::vector<shift_type>,
                                                       std::vector<shift_type>>& nodes,
                                       std::set<hash_type>&                      extras) {
    return std::make_pair(filter_nodes(graph, nodes.first, extras),
                          filter_nodes(graph, nodes.second, extras));
}


template <class GraphType>
typename Traverse<GraphType>::EndState
Traverse<GraphType>::dBG::traverse_left(GraphType *        graph,
                                        const std::string& seed,
                                        Path&              path,
                                        std::set<hash_type>&  mask) {

    this->set_cursor(seed);
    auto seed_hash = this->get();
    if (!graph->query(seed_hash)) {
        return {State::BAD_SEED, seed_hash};
    } 
    this->get_cursor(path);
    return traverse_left(graph, path, mask);
}


template <class GraphType>
typename Traverse<GraphType>::EndState
Traverse<GraphType>::dBG::traverse_left(GraphType * graph,
                                        Path&       path,
                                        std::set<hash_type>& mask) {

    hash_type prev_hash = this->get();
    this->seen.clear();
    this->seen.insert(prev_hash);

    // take the first step without check for reverse d-nodes
    std::pair<State, shift_type> step = step_left(graph);

    if (step.first != State::STEP_LEFT) {
        return {step.first, this->get()};
    } else {
        this->shift_left(step.second.symbol);
        path.push_front(step.second.symbol);
        this->seen.insert(step.second.hash);
    }

    while (1) {
        if (out_degree(graph) > 1) {
            pdebug("Stop: reverse d-node");
            path.pop_front();
            return {State::DECISION_RC, prev_hash};
        }

        prev_hash = this->get();
        std::pair<State, shift_type> step = step_left(graph);

        if (step.first != State::STEP_LEFT) {
            return {step.first, this->get()};
        } else {
            this->shift_left(step.second.symbol);
            path.push_front(step.second.symbol);
            this->seen.insert(step.second.hash);
        }
    }
}


template <class GraphType>
typename Traverse<GraphType>::EndState
Traverse<GraphType>::dBG::traverse_right(GraphType *        graph,
                                         const std::string& seed,
                                         Path&              path,
                                         std::set<hash_type>&  mask) {

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
                                         std::set<hash_type>& mask) {

    hash_type prev_hash = this->get();
    this->seen.clear();
    this->seen.insert(prev_hash);

    std::pair<State, shift_type> step = step_right(graph);

    if (step.first != State::STEP_RIGHT) {
        return {step.first, this->get()};
    } else {
        this->shift_right(step.second.symbol);
        path.push_back(step.second.symbol);
        this->seen.insert(step.second.hash);
    }
    
    while (1) {
        if (in_degree(graph) > 1) {
            path.pop_back();
            return {State::DECISION_RC, prev_hash};
        }

        // save the current cursor as prev_hash in case
        // we need to move backwards because of a reverse decision
        // k-mer
        prev_hash = this->get();
        std::pair<State, shift_type> step = step_right(graph);

        if (step.first != State::STEP_RIGHT) {
            return {step.first, this->get()};
        } else {
            this->shift_right(step.second.symbol);
            path.push_back(step.second.symbol);
            this->seen.insert(step.second.hash);
        }
    }
}


template <class GraphType>
std::pair<typename Traverse<GraphType>::EndState, typename Traverse<GraphType>::EndState>
Traverse<GraphType>::dBG::traverse(GraphType *        graph,
                                   const std::string& seed,
                                   Path&              path,
                                   std::set<hash_type>   mask) {

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

template struct boink::Traverse<boink::dBG<boink::storage::SparseppSetStorage,
                                           boink::hashing::UKHS::LazyShifter>>;
template struct boink::Traverse<boink::dBG<boink::storage::BitStorage,
                                           boink::hashing::UKHS::LazyShifter>>;
template struct boink::Traverse<boink::dBG<boink::storage::ByteStorage,
                                           boink::hashing::UKHS::LazyShifter>>;
template struct boink::Traverse<boink::dBG<boink::storage::NibbleStorage,
                                           boink::hashing::UKHS::LazyShifter>>;
template struct boink::Traverse<boink::dBG<boink::storage::QFStorage,
                                           boink::hashing::UKHS::LazyShifter>>;


#undef pdebug
