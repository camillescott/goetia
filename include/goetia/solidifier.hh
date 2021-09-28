/**
 * (c) Camille Scott, 2019
 * File   : solidifier.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 12.03.2020
 */

#ifndef GOETIA_SOLIDIFIER_HH
#define GOETIA_SOLIDIFIER_HH

#include "goetia/processors.hh"
#include "goetia/dbg.hh"
#include "goetia/storage/storage.hh"


namespace goetia {

template <class T>
struct StreamingSolidFilter;

template <template <class, class> class GraphType, class StorageType, class ShifterType>
struct StreamingSolidFilter<GraphType<StorageType, ShifterType>> {

    typedef GraphType<StorageType, ShifterType> graph_type;
    typedef ShifterType                         shifter_type;

    typedef typename shifter_type::alphabet     alphabet;
    typedef typename shifter_type::hash_type    hash_type;
   	typedef typename hash_type::value_type      value_type;
    typedef typename shifter_type::kmer_type    kmer_type;

    class Filter {
      private:
        std::shared_ptr<graph_type> dbg;

      public:
        const uint16_t    K;
        const float       min_prop_solid;
        const uint32_t    solid_threshold;

        Filter(std::shared_ptr<graph_type> dbg,
               const float                 min_prop_solid=0.75,
               const uint32_t              solid_threshold=1)
            : dbg(dbg),
              K(dbg->K),
              min_prop_solid(min_prop_solid),
              solid_threshold(solid_threshold)
        {
        }

        static std::shared_ptr<Filter> build(std::shared_ptr<graph_type> dbg,
                                             const float                 min_prop_solid=0.75,
                                             const uint32_t              solid_threshold=1) {
            return std::make_shared<Filter>(dbg, min_prop_solid, solid_threshold);
        }

        std::tuple<bool, uint64_t> filter_sequence(const std::string& sequence) {
            auto counts = dbg->insert_and_query_sequence(sequence);
            uint32_t n_not_solid = 0;
            uint32_t n_kmers = sequence.length() - K + 1;

            for (const auto& count : counts) {
                if (count < solid_threshold) {
                    ++n_not_solid;
                }
            }

            if (((float)n_not_solid / (float)n_kmers) > (1.0 - min_prop_solid)) {
                return {false, n_kmers};
            }
            return {true, n_kmers};
        }

    };

    using Processor = FilterProcessor<Filter>;

};

}

extern template class goetia::StreamingSolidFilter<goetia::dBG<goetia::storage::BitStorage, goetia::hashing::FwdLemireShifter>>;
extern template class goetia::StreamingSolidFilter<goetia::dBG<goetia::storage::QFStorage, goetia::hashing::FwdLemireShifter>>;


#endif
