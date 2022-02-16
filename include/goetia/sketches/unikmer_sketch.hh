/**
 * (c) Camille Scott, 2019
 * File   : ukhs_sketch.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */

#ifndef GOETIA_UKHS_SIGNATURE_HH
#define GOETIA_UKHS_SIGNATURE_HH

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <memory>
#include <set>
#include <vector>

#include "goetia/goetia.hh"
#include "goetia/processors.hh"

#include "goetia/hashing/kmeriterator.hh"
#include "goetia/hashing/hashextender.hh"
#include "goetia/hashing/ukhs.hh"
#include "goetia/hashing/unikmershifter.hh"
#include "goetia/hashing/canonical.hh"
//#include "goetia/sketches/hllcounter.hh"
#include "goetia/storage/storage_types.hh"

#include "goetia/pdbg.hh"

namespace goetia {


class IncompatibleSketch : public GoetiaException {
public:
    explicit IncompatibleSketch(const std::string& msg = "Incompatible sketches.")
        : GoetiaException(msg) { }
};


template <class StorageType, class HashType>
struct UnikmerSketch {

    typedef StorageType                                 storage_type;
    typedef StorageTraits<StorageType>                  storage_traits;
    typedef UnikmerLemirePolicy<HashType>      shift_policy;
    typedef HashShifter<shift_policy>          shifter_type;
    typedef typename shifter_type::ukhs_type            ukhs_type;
    typedef PdBG<StorageType, shifter_type>             pdbg_type;

    _goetia_model_typedefs_from_shiftertype(shifter_type)

    class Sketch {

    protected:

        std::shared_ptr<ukhs_type> ukhs_map;
        std::shared_ptr<pdbg_type> sketch;

    public:

        const uint16_t W;
        const uint16_t K;

        explicit Sketch(uint16_t W,
                           uint16_t K,
                           std::shared_ptr<ukhs_type> ukhs_map)
            : Sketch(W, K, ukhs_map, storage_traits::default_params)
        {
        }

        explicit Sketch(uint16_t W,
                           uint16_t K,
                           std::shared_ptr<ukhs_type> ukhs_map,
                           const typename storage_traits::params_type& storage_params)
            : W         (W),
              K         (K),
              ukhs_map  (ukhs_map)
        {
            sketch = std::make_shared<pdbg_type>(W,
                                                    K,
                                                    ukhs_map,
                                                    storage_params);
        }

        static std::shared_ptr<Sketch> build(uint16_t W,
                                                uint16_t K,
                                                std::shared_ptr<ukhs_type> ukhs_map) {
            return std::make_shared<Sketch>(W, K, ukhs_map, storage_traits::default_params);
        }

        static std::shared_ptr<Sketch> build(uint16_t W,
                                                uint16_t K,
                                                std::shared_ptr<ukhs_type> ukhs_map,
                                                const typename storage_traits::params_type& storage_params) {
            return std::make_shared<Sketch>(W, K, ukhs_map, storage_params);
        }

        inline void insert(const std::string& kmer) {
            sketch->insert(kmer);
        }

        inline size_t insert_sequence(const std::string& sequence) {
            sketch->insert_sequence(sequence);
            return sequence.length() - K + 1;
        }

        size_t get_size() const {
            return sketch->n_partitions();
        }

        std::vector<size_t> get_sketch_as_vector() {
            return sketch->get_partition_counts();
        }

        void * get_sketch_as_buffer() {
            return sketch->get_partition_counts_as_buffer();
        }

        uint64_t get_n_kmers() const {
            size_t n_kmers = 0;
            for (const auto& c : sketch->get_partition_counts()) {
                n_kmers += c;
            }
            return n_kmers;
        }

        /*
        std::set<hash_type> intersection(Sketch * other) {
            if (other->get_size() != this->get_size() or
                other->bucket_K   != this->bucket_K or
                other->K()        != this->K()) {
                
                throw IncompatibleSketch("Error: Sketchs not compatible");
            }
        }
        */
    };

    using Processor = InserterProcessor<Sketch>; 
 

};

extern template class goetia::UnikmerSketch<goetia::BitStorage, goetia::Hash<uint64_t>>;
extern template class goetia::UnikmerSketch<goetia::BitStorage, goetia::Canonical<uint64_t>>;

extern template class goetia::UnikmerSketch<goetia::SparseppSetStorage, goetia::Hash<uint64_t>>;
extern template class goetia::UnikmerSketch<goetia::SparseppSetStorage, goetia::Canonical<uint64_t>>;

extern template class goetia::UnikmerSketch<goetia::PHMapStorage, goetia::Hash<uint64_t>>;
extern template class goetia::UnikmerSketch<goetia::PHMapStorage, goetia::Canonical<uint64_t>>;

extern template class goetia::UnikmerSketch<goetia::BTreeStorage, goetia::Hash<uint64_t>>;
extern template class goetia::UnikmerSketch<goetia::BTreeStorage, goetia::Canonical<uint64_t>>;

extern template class goetia::UnikmerSketch<goetia::ByteStorage, goetia::Hash<uint64_t>>;
extern template class goetia::UnikmerSketch<goetia::ByteStorage, goetia::Canonical<uint64_t>>;

extern template class goetia::UnikmerSketch<goetia::NibbleStorage, goetia::Hash<uint64_t>>;
extern template class goetia::UnikmerSketch<goetia::NibbleStorage, goetia::Canonical<uint64_t>>;

extern template class goetia::UnikmerSketch<goetia::QFStorage, goetia::Hash<uint64_t>>;
extern template class goetia::UnikmerSketch<goetia::QFStorage, goetia::Canonical<uint64_t>>;

//extern template class goetia::UnikmerSketch<goetia::HLLStorage, goetia::Hash<uint64_t>>;
//extern template class goetia::UnikmerSketch<goetia::HLLStorage, goetia::Canonical<uint64_t>>;

}


#endif
