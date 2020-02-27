/**
 * (c) Camille Scott, 2019
 * File   : kmeriterator.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 23.07.2019
 */
/* kmeriterator.hh -- sequence k-mer iterators
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef GOETIA_KMERITERATOR_HH
#define GOETIA_KMERITERATOR_HH

#include <string>

#include "goetia/goetia.hh"
#include "goetia/hashing/hashshifter.hh"
#include "goetia/sequences/exceptions.hh"

#include "goetia/hashing/shifter_types.hh"


namespace goetia {
namespace hashing {

template <class ShifterType>
class KmerIterator {

    const std::string _seq;
    unsigned int index;
    unsigned int length;
    bool _initialized, _shifter_owner;

public:

    typedef ShifterType                     shifter_type;
    typedef typename ShifterType::hash_type hash_type;
    const uint16_t                          K;

    ShifterType * shifter;

    template<typename... Args>
    __attribute__((visibility("default")))
    explicit KmerIterator(const std::string& seq, uint16_t K, Args&&... args)
        : K(K), 
          _seq(seq), 
          index(0), 
          _initialized(false), 
          _shifter_owner(true)
    {

        if (_seq.length() < K) {
            throw SequenceLengthException("Sequence must have length >= K");
        }
        shifter = new ShifterType(seq, K, std::forward<Args>(args)...);
    }

    __attribute__((visibility("default")))
    explicit KmerIterator(const std::string& seq, ShifterType * shifter)
        : K(shifter->K), 
          _seq(seq),
          index(0), 
          _initialized(false),
          _shifter_owner(false), 
          shifter(shifter) 
    {
        if (_seq.length() < K) {
            throw SequenceLengthException("Sequence must have length >= K");
        }
        shifter->hash_base(_seq);
    }

    __attribute__((visibility("default")))
    explicit KmerIterator(const std::string& seq, ShifterType& shifter_proto)
        : K(shifter_proto.K), 
          _seq(seq), 
          index(0), 
          _initialized(false), 
          _shifter_owner(true)
    {

        if (_seq.length() < K) {
            throw SequenceLengthException("Sequence must have length >= K");
        }
        shifter = new ShifterType(shifter_proto);
        shifter->hash_base(_seq);
    }

    __attribute__((visibility("default")))
    ~KmerIterator() {
        if (_shifter_owner) {
            delete shifter;
        }
    }

    __attribute__((visibility("default")))
    hash_type first()  {
        _initialized = true;
        index += 1;
        return shifter->get();
    }

    __attribute__((visibility("default")))
    hash_type next() {
        if (!_initialized) {
            return first();
        }

        if (done()) {
            throw InvalidCharacterException("past end of iterator");
        }

        auto ret = shifter->shift_right(_seq[index - 1], _seq[index + K - 1]);
        index += 1;

        return ret;
    }

    __attribute__((visibility("default")))
    bool done() const  {
        return (index + K > _seq.length());
    }

    __attribute__((visibility("default")))
    unsigned int get_start_pos() const {
        if (!_initialized) { return 0; }
        return index - 1;
    }

    __attribute__((visibility("default")))
    unsigned int get_end_pos() const {
        if (!_initialized) { return K; }
        return index + K - 1;
    }
};


extern template class KmerIterator<FwdLemireShifter>;
extern template class KmerIterator<CanLemireShifter>;

extern template class KmerIterator<FwdUnikmerShifter>;
extern template class KmerIterator<CanUnikmerShifter>;


}
}

#undef pdebug
#endif
