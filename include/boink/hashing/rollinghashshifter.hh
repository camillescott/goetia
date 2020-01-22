/**
 * (c) Camille Scott, 2019
 * File   : rollinghashshifter.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 06.08.2019
 */
/* rollinghashshifter.hh -- k-mer hash functions
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_ROLLINGHASHSHIFTER_HH
#define BOINK_ROLLINGHASHSHIFTER_HH

#include "boink/boink.hh"
#include "boink/meta.hh"
#include "boink/sequences/alphabets.hh"
#include "boink/hashing/canonical.hh"
#include "boink/hashing/hashshifter.hh"

#include "boink/hashing/rollinghash/cyclichash.h"

namespace boink::hashing {


template<class HashType>
struct RollingHashShifterBase {
protected:
    typedef typename HashType::value_type value_type;

    CyclicHash<value_type> hasher;

    RollingHashShifterBase(uint16_t K)
        : hasher(K) {}
};

template<>
struct RollingHashShifterBase<CanonicalModel<uint64_t>> {
protected:
    typedef typename CanonicalModel<uint64_t>::value_type value_type;

    CyclicHash<value_type> hasher;
    CyclicHash<value_type> rc_hasher;

    RollingHashShifterBase(uint16_t K)
        : hasher(K),
          rc_hasher(K) {}
};


template<class HashType = HashModel<uint64_t>>
class RollingHashShifter : public HashShifter<RollingHashShifter<HashType>,
                                              HashType,
                                              DNA_SIMPLE>,
                           public RollingHashShifterBase<HashType> {

    typedef HashShifter<RollingHashShifter<HashType>,
                        HashType,
                        DNA_SIMPLE> BaseShifter;
    typedef Tagged<RollingHashShifter<HashType>> tagged_type;

public:

    friend BaseShifter;

    using BaseShifter::NAME;
    using BaseShifter::OBJECT_ABI_VERSION;

    typedef typename BaseShifter::value_type value_type;
    typedef typename BaseShifter::hash_type  hash_type;
    typedef typename BaseShifter::kmer_type  kmer_type;
    typedef typename BaseShifter::alphabet   alphabet;


    RollingHashShifter(const std::string& start,
                       uint16_t K)
        : BaseShifter(start, K),
          RollingHashShifterBase<HashType>(K)
    {    
    }

    RollingHashShifter(uint16_t K)
        : BaseShifter(K),
          RollingHashShifterBase<HashType>(K)
    {
    }

    RollingHashShifter(RollingHashShifter const& other)
        : RollingHashShifter(other.K())
    {
    }

    __attribute__((visibility("default")))
    inline hash_type _hash_base(const char * sequence) {
        this->hasher.reset();
        for (uint16_t i = 0; i < this->_K; ++i) {
            if (sequence[i] == '\0') {
                throw SequenceLengthException("Encountered null terminator in k-mer!");
            }
            this->hasher.eat(sequence[i]);
        }
        return _get();
    }

    template<class It> __attribute__((visibility("default")))
    inline hash_type _hash_base(It begin, It end) {
        this->hasher.reset();
        while (begin != end) {
            this->hasher.eat(*begin);
            begin = std::next(begin);
        }

        return _get();
    }

    hash_type _get() {
        return {this->hasher.hashvalue};
    }

    static hash_type _hash(const char * sequence, const uint16_t K) {
        CyclicHash<value_type> tmp_hasher(K);
        for (uint16_t i = 0; i < K; ++i) {
            tmp_hasher.eat(sequence[i]);
        }
        return {tmp_hasher.hashvalue};
    }

    hash_type _shift_left(const char& in, const char& out) {
        this->hasher.reverse_update(in, out);
        return _get();
    }


    hash_type _shift_right(const char& out, const char& in) {
        this->hasher.update(out, in);
        return _get();
    }

};


//
// Canonical specialization decls (sigh)
//

template<>
inline RollingHashShifter<CanonicalModel<uint64_t>>
::RollingHashShifter(const std::string& start,
                     uint16_t K)
    : BaseShifter(start, K),
      RollingHashShifterBase(K)
{
}


template<>
inline RollingHashShifter<CanonicalModel<uint64_t>>
::RollingHashShifter(uint16_t K)
    : BaseShifter(K),
      RollingHashShifterBase(K)
{
}


template<>
inline CanonicalModel<uint64_t>
RollingHashShifter<CanonicalModel<uint64_t>>
::_get() {
    return {hasher.hashvalue, rc_hasher.hashvalue};
}


// static
template<>
inline CanonicalModel<uint64_t>
RollingHashShifter<CanonicalModel<uint64_t>>
::_hash(const char * sequence, const uint16_t K)  {
    CyclicHash<value_type> fw_hasher(K);
    CyclicHash<value_type> rc_hasher(K);

    for (uint16_t i = 0; i < K; ++i) {
        fw_hasher.eat(sequence[i]);
        rc_hasher.eat(alphabet::complement(sequence[K - i - 1]));
    }

    return {fw_hasher.hashvalue, rc_hasher.hashvalue};
}


template<>
inline CanonicalModel<uint64_t>
RollingHashShifter<CanonicalModel<uint64_t>>
::_hash_base(const char * sequence) {
    hasher.reset();
    rc_hasher.reset();

    for (uint16_t i = 0; i < this->_K; ++i) {
        if (sequence[i] == '\0') {
            throw SequenceLengthException("Encountered null terminator in k-mer!");
        }
        hasher.eat(sequence[i]);
        rc_hasher.eat(alphabet::complement(sequence[this->_K - i - 1]));
    }

    return _get();
}


template<>
template<class It>
inline CanonicalModel<uint64_t>
RollingHashShifter<CanonicalModel<uint64_t>>
::_hash_base(It begin, It end) {

    hasher.reset();
    rc_hasher.reset();
    --end; // the end() iter points to last position + 1

    size_t i = 0;
    while (i < this->K()) {
        hasher.eat(*(begin));
        rc_hasher.eat(alphabet::complement(*(end)));

        //std::cout << "eat: " << *begin << " fwd, " << *end << " rc ("
        //          << alphabet::complement(*(end)) << ")" << std::endl;

        ++begin;
        --end;

        ++i;
    }

    return _get();
}


template<>
inline CanonicalModel<uint64_t>
RollingHashShifter<CanonicalModel<uint64_t>>
::_shift_right(const char& out, const char& in) {
    hasher.update(out, in);
    rc_hasher.reverse_update(alphabet::complement(in),
                             alphabet::complement(out));
    return _get();
}


template<>
inline CanonicalModel<uint64_t>
RollingHashShifter<CanonicalModel<uint64_t>>
::_shift_left(const char& in, const char& out) {
    hasher.reverse_update(in, out);
    rc_hasher.update(alphabet::complement(out),
                     alphabet::complement(in));
    return _get();
}

typedef RollingHashShifter<HashModel<uint64_t>> FwdRollingShifter;
typedef RollingHashShifter<CanonicalModel<uint64_t>> CanRollingShifter;

extern template class RollingHashShifter<HashModel<uint64_t>>;
extern template class RollingHashShifter<CanonicalModel<uint64_t>>;

} 

#endif
