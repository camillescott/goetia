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

template<class GraphType>
template<Direction_t D>
size_t dBGWalker<GraphType>
::count_nodes(GraphType *                       graph,
              const std::vector<shift_type<D>>& extensions) {

    uint8_t n_found = 0;
    for (const auto& ext : extensions) {
        if(graph->query(ext.value())) {
            ++n_found;
        }
    }
    return n_found;
}


template<class GraphType>
template<Direction_t D>
size_t dBGWalker<GraphType>
::count_nodes(GraphType *                graph,
              const std::vector<shift_type<D>>& extensions,
              std::set<hash_type>&              extras) {

    uint8_t n_found = 0;
    for (const auto& ext : extensions) {
        if(graph->query(ext.value()) ||
           extras.count(ext.value())) {
            ++n_found;
        }
    }
    return n_found;
}


template<class GraphType>
template<Direction_t D>
auto dBGWalker<GraphType>
::reduce_nodes(GraphType *                       graph,
               const std::vector<shift_type<D>>& extensions)
-> shift_type<D> {

    uint8_t n_found = 0;
    shift_type<D> result;
    for (const auto& ext : extensions) {
        if(graph->query(ext.value())) {
            ++n_found;
            if (n_found > 1) {
                return n_found;
            }
            result = ext;
        }
    }
    return result;
}


template<class GraphType>
template<Direction_t D>
auto dBGWalker<GraphType>
::reduce_nodes(GraphType *                       graph,
               const std::vector<shift_type<D>>& extensions,
               std::set<hash_type>&              extra)
-> shift_type<D> {
    uint8_t n_found = 0;
    shift_type<D> result;
    for (const auto& ext : extensions ) {
        if(graph->query(ext.value()) ||
           extra.count(ext.value())) {
            ++n_found;
            if (n_found > 1) {
                return n_found;
            }
            result = ext;
        }
    }
    return result;
}


template<class GraphType>
template<Direction_t D>
auto dBGWalker<GraphType>
::filter_nodes(GraphType * graph,
               const std::vector<shift_type<D>>& extensions)
-> std::vector<shift_type<D>> {

    std::vector<shift_type<D>> result;
    for (const auto& ext : extensions) {
        if (graph->query(ext.value())) {
            result.push_back(ext);
        }
    }
    return result;
}


template<class GraphType>
template<Direction_t D>
auto dBGWalker<GraphType>
::filter_nodes(GraphType *                       graph,
               const std::vector<shift_type<D>>& extensions,
               std::set<hash_type>&              extra)
-> std::vector<shift_type<D>> {

    std::vector<shift_type<D>> result;
    for (const auto& ext : extensions) {
        if (graph->query(ext.value()) ||
            extra.count(ext.value())) {
            result.push_back(ext);
        }
    }
    return result;
}


template <class GraphType>
auto dBGWalker<GraphType>::filter_nodes(GraphType *            graph,
                                        const shift_pair_type& extensions)
-> shift_pair_type {
    return std::make_pair(filter_nodes(graph, extensions.first),
                          filter_nodes(graph, extensions.second));
}


template <class GraphType>
auto dBGWalker<GraphType>
::filter_nodes(GraphType *            graph,
               const shift_pair_type& extensions,
               std::set<hash_type>&   extras)
-> shift_pair_type {
    return std::make_pair(filter_nodes(graph, extensions.first, extras),
                          filter_nodes(graph, extensions.second, extras));
}


template <class GraphType>
auto dBGWalker<GraphType>
::walk_left(GraphType *          graph,
            const std::string&   seed,
            std::set<hash_type>& mask)
-> Walk<hashing::DIR_LEFT> {

    this->set_cursor(seed);
    auto seed_hash = this->get();
    if (!graph->query(seed_hash)) {
        Walk<hashing::DIR_LEFT> walk{kmer_type{seed_hash, this->get_cursor()},
                                     {},
                                     State::BAD_SEED};
        return walk;
    } 
    return walk_left(graph, mask);
}


template <class GraphType>
auto dBGWalker<GraphType>
::walk_left(GraphType *          graph,
            std::set<hash_type>& mask)
-> Walk<hashing::DIR_LEFT> {

    Walk<hashing::DIR_LEFT> walk;

    hash_type start_hash = this->get();
    this->seen.clear();
    this->seen.insert(start_hash);
    walk.start.hash = start_hash;
    walk.start.kmer = this->get_cursor();

    // take the first step without check for reverse d-nodes
    auto step = step_left(graph, mask);

    if (step.first != State::STEP) {
        walk.end_state = step.first;
        return walk;
    } else {
        walk.path.push_back(step.second.front());
    }

    while (1) {
        if (out_degree(graph) > 1) {
            pdebug("Stop: reverse d-node");
            walk.path.pop_back();
            walk.end_state = State::DECISION_RC;
            this->seen.erase(this->get());
            return std::move(walk);
        }

        step = step_left(graph, mask);

        if (step.first != State::STEP) {
            walk.end_state = step.first;
            
            return std::move(walk);
        } else {
            walk.path.push_back(step.second.front());
        }
    }
}


template <class GraphType>
auto dBGWalker<GraphType>
::walk_right(GraphType *            graph,
             const std::string&     seed,
             std::set<hash_type>&  mask)
-> Walk<hashing::DIR_RIGHT> {

    this->set_cursor(seed);
    auto seed_hash = this->get();
    if (!graph->query(seed_hash)) {
        Walk<hashing::DIR_RIGHT> walk{kmer_type{seed_hash, this->get_cursor()},
                                       {},
                                       State::BAD_SEED};
        return walk;
    }
    return walk_right(graph, mask);
}


template <class GraphType>
auto dBGWalker<GraphType>
::walk_right(GraphType *           graph,
             std::set<hash_type>& mask)
-> Walk<hashing::DIR_RIGHT> {

    Walk<hashing::DIR_RIGHT> walk;

    hash_type start_hash = this->get();
    this->seen.clear();
    this->seen.insert(start_hash);
    walk.start.hash = start_hash;
    walk.start.kmer = this->get_cursor();

    auto step = step_right(graph, mask);

    if (step.first != State::STEP) {
        walk.end_state = step.first;
        return walk;
    } else {
        walk.path.push_back(step.second.front());
    }
    
    while (1) {
        if (in_degree(graph) > 1) {
            walk.path.pop_back();
            walk.end_state = State::DECISION_RC;
            this->seen.erase(this->get());
            return std::move(walk);
        }

        step = step_right(graph, mask);

        if (step.first != State::STEP) {
            walk.end_state = step.first;
            return std::move(walk);
        } else {
            walk.path.push_back(step.second.front());
        }
    }
}


template <class GraphType>
auto  dBGWalker<GraphType>
::walk(GraphType *          graph,
       const std::string&   seed,
       std::set<hash_type> mask) 
-> walk_pair_type {

    this->set_cursor(seed);
    auto seed_hash = this->get();
    if (!graph->query(seed_hash)) {
        kmer_type start{seed_hash, this->get_cursor()};
        return {{start, {}, State::BAD_SEED},
                {start, {}, State::BAD_SEED}};
    }

    auto left_walk = walk_left(graph, mask);
    this->set_cursor(seed);
    auto right_walk = walk_right(graph, mask);
    return {left_walk, right_walk};
}

template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>>;
template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::CanRollingShifter>>;
template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>>;
template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>>;


template class dBGWalker<dBG<storage::BitStorage, hashing::FwdRollingShifter>>;
template class dBGWalker<dBG<storage::BitStorage, hashing::CanRollingShifter>>;
template class dBGWalker<dBG<storage::BitStorage, hashing::FwdUnikmerShifter>>;
template class dBGWalker<dBG<storage::BitStorage, hashing::CanUnikmerShifter>>;

template class dBGWalker<dBG<storage::ByteStorage, hashing::FwdRollingShifter>>;
template class dBGWalker<dBG<storage::ByteStorage, hashing::CanRollingShifter>>;
template class dBGWalker<dBG<storage::ByteStorage, hashing::FwdUnikmerShifter>>;
template class dBGWalker<dBG<storage::ByteStorage, hashing::CanUnikmerShifter>>;

template class dBGWalker<dBG<storage::NibbleStorage, hashing::FwdRollingShifter>>;
template class dBGWalker<dBG<storage::NibbleStorage, hashing::CanRollingShifter>>;
template class dBGWalker<dBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>>;
template class dBGWalker<dBG<storage::NibbleStorage, hashing::CanUnikmerShifter>>;

template class dBGWalker<dBG<storage::QFStorage, hashing::FwdRollingShifter>>;
template class dBGWalker<dBG<storage::QFStorage, hashing::CanRollingShifter>>;
template class dBGWalker<dBG<storage::QFStorage, hashing::FwdUnikmerShifter>>;
template class dBGWalker<dBG<storage::QFStorage, hashing::CanUnikmerShifter>>;

} // namespace boink


#undef pdebug
