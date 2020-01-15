/* rollinghashshifter.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "boink/hashing/rollinghashshifter.hh"

namespace boink::hashing {

template<>
RollingHashShifter<CanonicalModel<uint64_t>>
::RollingHashShifter(const std::string& start,
                     uint16_t K)
    : BaseShifter(start, K),
      RollingHashShifterBase(K)
{
}


template<>
RollingHashShifter<CanonicalModel<uint64_t>>
::RollingHashShifter(uint16_t K)
    : BaseShifter(K),
      RollingHashShifterBase(K)
{
}


template<>
CanonicalModel<uint64_t>
RollingHashShifter<CanonicalModel<uint64_t>>
::_get() {
    return {hasher.hashvalue, rc_hasher.hashvalue};
}


// static
template<>
CanonicalModel<uint64_t>
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
CanonicalModel<uint64_t>
RollingHashShifter<CanonicalModel<uint64_t>>
::_hash_base(const char * sequence) {
    hasher.reset();
    rc_hasher.reset();

    for (uint16_t i = 0; i < this->_K; ++i) {
        hasher.eat(sequence[i]);
        rc_hasher.eat(alphabet::complement(sequence[this->_K - i - 1]));
    }

    return _get();
}


template<>
CanonicalModel<uint64_t>
RollingHashShifter<CanonicalModel<uint64_t>>
::_shift_right(const char& in, const char& out) {
    hasher.update(out, in);
    rc_hasher.reverse_update(alphabet::complement(in),
                             alphabet::complement(out));
    return _get();
}


template<>
CanonicalModel<uint64_t>
RollingHashShifter<CanonicalModel<uint64_t>>
::_shift_left(const char& in, const char& out) {
    hasher.reverse_update(in, out);
    rc_hasher.update(alphabet::complement(out),
                     alphabet::complement(in));
    return _get();
}

template class RollingHashShifter<HashModel<uint64_t>>;
template class RollingHashShifter<CanonicalModel<uint64_t>>;

}

namespace boink {

template <> const std::string
    Tagged<hashing::RollingHashShifter<hashing::HashModel<uint64_t>>>
    ::NAME = "Boink::CyclicHash<uint64_t>_NonRandom_Fwd";

template <> const std::string 
    Tagged<hashing::RollingHashShifter<hashing::CanonicalModel<uint64_t>>>
    ::NAME = "Boink::CyclicHash<uint64_t>_NonRandom_Can";

}
