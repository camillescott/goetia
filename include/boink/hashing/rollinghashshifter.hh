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

#include "boink/hashing/rollinghash/cyclichash.h"

namespace boink::hashing {


template<class HashType>
struct RollingHashShifterBase {
protected:
    typedef typename HashType::value_type value_type;

    CyclicHash<value_type> hasher;

    explicit RollingHashShifterBase(uint16_t K)
        : hasher(K) {
        
        //std::cout << "END RollingBase(K) ctor " << this << std::endl;
    }

    RollingHashShifterBase() = delete;
};

template<>
struct RollingHashShifterBase<Canonical<uint64_t>> {
protected:
    typedef typename Canonical<uint64_t>::value_type value_type;

    CyclicHash<value_type> hasher;
    CyclicHash<value_type> rc_hasher;

    explicit RollingHashShifterBase(uint16_t K)
        : hasher(K),
          rc_hasher(K) {}
};


template<typename HashType,
         typename Alphabet = DNA_SIMPLE>
class RollingHashShifter : public RollingHashShifterBase<HashType> {

    typedef Tagged<RollingHashShifter<HashType, Alphabet>> tagged_type;

public:

    typedef HashType                       hash_type;
    typedef typename hash_type::value_type value_type;
    typedef Kmer<hash_type>           kmer_type;
    typedef Alphabet                       alphabet;
    static constexpr bool has_kmer_span = false;

    const uint16_t K;

    __attribute__((visibility("default")))
    inline hash_type hash_base_impl(const char * sequence) {
        this->hasher.reset();
        for (uint16_t i = 0; i < K; ++i) {
            if (sequence[i] == '\0') {
                throw SequenceLengthException("Encountered null terminator in k-mer!");
            }
            this->hasher.eat(sequence[i]);
        }
        return get_impl();
    }

    template<class It> __attribute__((visibility("default")))
    inline hash_type hash_base_impl(It begin, It end) {
        this->hasher.reset();
        while (begin != end) {
            this->hasher.eat(*begin);
            begin = std::next(begin);
        }

        return get_impl();
    }

    hash_type get_impl() {
        //std::cout << "RollingHash: get_impl " << this->hasher.hashvalue << std::endl;;
        return {this->hasher.hashvalue};
    }

    hash_type shift_left_impl(const char& in, const char& out) {
        this->hasher.reverse_update(in, out);
        return get_impl();
    }


    hash_type shift_right_impl(const char& out, const char& in) {
        this->hasher.update(out, in);
        return get_impl();
    }

protected:

    explicit RollingHashShifter(uint16_t K)
        : RollingHashShifterBase<HashType>(K),
          K(K)
    {   
        //std::cout << "END RollingHashShifter(K) ctor " << this << " / " << static_cast<RollingHashShifterBase<HashType>*>(this) << std::endl;
    }

    explicit RollingHashShifter(const RollingHashShifter& other)
        : RollingHashShifter(other.K)
    {
        //std::cout << "END RollingHashShifter(other) ctor " << this << " / " << static_cast<RollingHashShifterBase<HashType>*>(this) << std::endl;
    }

    RollingHashShifter() = delete;

};


template<>
inline RollingHashShifter<Canonical<uint64_t>>
::RollingHashShifter(uint16_t K)
    : RollingHashShifterBase<Canonical<uint64_t>>(K),
      K(K)
{
}

template<>
inline RollingHashShifter<Canonical<uint64_t>>
::RollingHashShifter(const RollingHashShifter& other)
    : RollingHashShifterBase<Canonical<uint64_t>>(other.K),
      K(other.K)
{
}


template<>
inline Canonical<uint64_t>
RollingHashShifter<Canonical<uint64_t>>
::get_impl() {
    return {hasher.hashvalue, rc_hasher.hashvalue};
}

template<>
template<class It>
inline Canonical<uint64_t>
RollingHashShifter<Canonical<uint64_t>>
::hash_base_impl(It begin, It end) {

    hasher.reset();
    rc_hasher.reset();
    --end; // the end() iter points to last position + 1

    size_t i = 0;
    while (i < K) {
        hasher.eat(*(begin));
        rc_hasher.eat(alphabet::complement(*(end)));

        //std::cout << "eat: " << *begin << " fwd, " << *end << " rc ("
        //          << alphabet::complement(*(end)) << ")" << std::endl;

        ++begin;
        --end;

        ++i;
    }

    return get_impl();
}


template<>
inline Canonical<uint64_t>
RollingHashShifter<Canonical<uint64_t>>
::hash_base_impl(const char * sequence) {
    hasher.reset();
    rc_hasher.reset();

    for (uint16_t i = 0; i < K; ++i) {
        if (sequence[i] == '\0') {
            throw SequenceLengthException("Encountered null terminator in k-mer!");
        }
        hasher.eat(sequence[i]);
        rc_hasher.eat(alphabet::complement(sequence[K - i - 1]));
    }

    return get_impl();
}


template<>
inline Canonical<uint64_t>
RollingHashShifter<Canonical<uint64_t>>
::shift_right_impl(const char& out, const char& in) {
    hasher.update(out, in);
    rc_hasher.reverse_update(alphabet::complement(in),
                             alphabet::complement(out));
    return get_impl();
}


template<>
inline Canonical<uint64_t>
RollingHashShifter<Canonical<uint64_t>>
::shift_left_impl(const char& in, const char& out) {
    hasher.reverse_update(in, out);
    rc_hasher.update(alphabet::complement(out),
                     alphabet::complement(in));
    return get_impl();
}

typedef RollingHashShifter<Hash<uint64_t>> FwdLemirePolicy;
typedef RollingHashShifter<Canonical<uint64_t>> CanLemirePolicy;

extern template class RollingHashShifter<Hash<uint64_t>>;
extern template class RollingHashShifter<Canonical<uint64_t>>;

} 

#endif
