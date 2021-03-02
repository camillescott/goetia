/**
 * (c) Camille Scott, 2019
 * File   : cdbg.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 03.09.2019
 */

#include "goetia/cdbg/cdbg.hh"

#include "goetia/dbg.hh"
#include "goetia/hashing/rollinghashshifter.hh"
#include "goetia/storage/storage_types.hh"

namespace goetia {
namespace cdbg {

template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
auto
cDBG<GraphType<StorageType, ShifterType>>
::compute_connected_component_metrics(std::shared_ptr<cDBG<GraphType<StorageType, ShifterType>>::Graph> cdbg,
                                      size_t sample_size)
        -> std::tuple<size_t, size_t, size_t, std::vector<size_t>> {

        auto time_start = std::chrono::system_clock::now();

        metrics::ReservoirSample<size_t> component_size_sample(sample_size);
        size_t max_component = 0;
        size_t min_component = std::numeric_limits<size_t>::max();
        auto components = cdbg->find_connected_components();

        for (auto id_comp_pair : components) {
            size_t component_size = id_comp_pair.second.size();
            component_size_sample.sample(component_size);
            max_component = (component_size > max_component) ? component_size : max_component;
            min_component = (component_size < min_component) ? component_size : min_component;
        }

        auto time_elapsed = std::chrono::system_clock::now() - time_start;
        _cerr("Finished recomputing components. Elapsed time: " <<
              std::chrono::duration<double>(time_elapsed).count());

        return {components.size(), min_component, max_component, component_size_sample.get_result()};
    }

}

template class cdbg::cDBG<dBG<storage::BitStorage, hashing::FwdLemireShifter>>;
template class cdbg::cDBG<dBG<storage::BitStorage, hashing::CanLemireShifter>>;

template class cdbg::cDBG<dBG<storage::SparseppSetStorage, hashing::FwdLemireShifter>>;
template class cdbg::cDBG<dBG<storage::SparseppSetStorage, hashing::CanLemireShifter>>;

template class cdbg::cDBG<dBG<storage::ByteStorage, hashing::FwdLemireShifter>>;
template class cdbg::cDBG<dBG<storage::ByteStorage, hashing::CanLemireShifter>>;

template class cdbg::cDBG<dBG<storage::NibbleStorage, hashing::FwdLemireShifter>>;
template class cdbg::cDBG<dBG<storage::NibbleStorage, hashing::CanLemireShifter>>;

template class cdbg::cDBG<dBG<storage::QFStorage, hashing::FwdLemireShifter>>;
template class cdbg::cDBG<dBG<storage::QFStorage, hashing::CanLemireShifter>>;
}
