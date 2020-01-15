/**
 * (c) Camille Scott, 2019
 * File   : ukhs.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 07.08.2019
 */
/* ukhs.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include <limits>

#include "boink/hashing/ukhs.hh"


namespace boink::hashing {

template <class ShifterType>
UKHS<ShifterType>::UKHS(uint16_t W,
                        uint16_t K,
                        std::vector<std::string>& ukhs)
    : KmerClient  (K),
      _W          (W)
{
    if (ukhs.front().size() != K) {
        throw BoinkException("K does not match k-mer size from provided UKHS");
    }

    ShifterType hasher(K);
    uint64_t pid = 0;
    for (auto unikmer : ukhs) {
        hash_type h = hasher.hash(unikmer);

        if (!pmap.count(h.value())) {
            hashes.push_back(h);
            pmap[h.value()] = pid;
            ++pid;
        }
    }
}


template <class ShifterType>
std::optional<Partitioned<typename ShifterType::hash_type>>
UKHS<ShifterType>::query(hash_type unikmer_hash) {

    auto search = pmap.find(unikmer_hash.value());
    if (search != pmap.end()) {
        return {unikmer_hash, search->second};
    } else {
        return {};
    }
}

}

