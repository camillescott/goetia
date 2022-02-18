/**
 * (c) Camille Scott, 2021
 * File   : streamhasher.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 17.02.2022
 */

#ifndef GOETIA_STREAMHASHER_HH
#define GOETIA_STREAMHASHER_HH

#include <memory>

#include "goetia/hashing/canonical.hh"
#include "goetia/hashing/kmeriterator.hh"
#include "goetia/processors.hh"

namespace goetia {

template <class ShifterType>
struct StreamHasher {

    typedef ShifterType                         shifter_type;
    typedef typename shifter_type::alphabet     alphabet;
    typedef typename shifter_type::hash_type    hash_type;
   	typedef typename hash_type::value_type      value_type;
    typedef typename shifter_type::kmer_type    kmer_type;

    class Hasher {
      public:
        const uint16_t K;

        Hasher(uint16_t K)
            : K(K)
        {
        }

        static std::shared_ptr<Hasher> build(uint16_t K) {
            return std::make_shared<Hasher>(K);
        }

        uint64_t insert_sequence(const std::string& sequence) {
            KmerIterator<ShifterType> iter(sequence, K);
            while(!iter.done()) {
                hash_type h = iter.next();
            }
            return sequence.size() - K + 1;
        }
    };

    using Processor = InserterProcessor<Hasher>;
};

extern template class StreamHasher<FwdLemireShifter>;
extern template class StreamHasher<CanLemireShifter>;

}

#endif
