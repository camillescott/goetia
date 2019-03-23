/* ukhs.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "boink/hashing/ukhs.hh"

namespace boink {
namespace hashing {

template <>
template <>
typename hash_return<PartitionedHash>::type
KmerIterator<UKHShifter>::first<PartitionedHash>() {
    _initialized = true;
    index += 1;
    shifter->reset_unikmers();
    return std::make_pair<hash_t, uint64_t>(shifter->get(), shifter->get_back_partition());
}


template <>
template <>
typename hash_return<PartitionedHash>::type
KmerIterator<UKHShifter>::next<PartitionedHash>() {
    if (!_initialized) {
        return this->first<PartitionedHash>();
    }

    if (done()) {
        throw InvalidCharacterException("past end of iterator");
    }

    shifter->shift_right(_seq[index + _K - 1]);
    index += 1;

    return std::make_pair<hash_t, uint64_t>(shifter->get(), shifter->get_back_partition());
}

}
}
