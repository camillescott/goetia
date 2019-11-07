/* rollinghashshifter.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/kmeriterator.hh"

namespace boink {

template <> const std::string Tagged<hashing::RollingHashShifter>::NAME = "Boink::CyclicHash<uint64_t>::NonRandom";

}

template class boink::hashing::KmerIterator<boink::hashing::RollingHashShifter>;
