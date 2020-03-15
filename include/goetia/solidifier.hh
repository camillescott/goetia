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

        Filter(std::shared_ptr<graph_type> dbg,
                   const float             min_prop_solid=0.75)
            : dbg(dbg),
              K(dbg->K),
              min_prop_solid(min_prop_solid)
        {
        }

        static std::shared_ptr<Filter> build(std::shared_ptr<graph_type> dbg,
                                             const float min_prop_solid=0.75) {
            return std::make_shared<Filter>(dbg, min_prop_solid);
        }

        bool filter_sequence(const std::string& sequence) {
            uint64_t n_not_solid = dbg->insert_sequence(sequence);
            uint64_t n_kmers = sequence.length() - K + 1;
            
            if (((float)n_not_solid / (float)n_kmers) > (1.0 - min_prop_solid)) {
                return false;
            }
            return true;
        }

    };

    using Processor = FilterProcessor<Filter>;

};

extern template class StreamingSolidFilter<dBG<storage::BitStorage, hashing::FwdLemireShifter>>;
extern template class StreamingSolidFilter<dBG<storage::QFStorage, hashing::FwdLemireShifter>>;

}
#endif
