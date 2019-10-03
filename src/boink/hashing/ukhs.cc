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

#include "boink/hashing/ukhs.hh"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#include "bbhash/BooPHF.h"
#pragma GCC diagnostic pop


namespace boink {
namespace hashing {

UKHS::Map::Map(uint16_t W,
               uint16_t K,
               std::vector<std::string>& ukhs)
    : KmerClient  (K),
      ukhs_revmap (ukhs.size()),
      _W          (W)
{
    if (ukhs.front().size() != K) {
        throw BoinkException("K does not match k-mer size from provided UKHS");
    }

    //std::cerr << "Building MPHF with "
    //          << ukhs.size() << " k-mers."
    //          << std::endl;

    for (auto unikmer : ukhs) {
        ukhs_hashes.push_back(hash_cyclic(unikmer, K));
    }
    bphf = std::make_unique<boophf_t>(ukhs_hashes.size(),
                                      ukhs_hashes,
                                      1,
                                      2.0,
                                      false,
                                      false);
    for (auto unikmer_hash : ukhs_hashes) {
        ukhs_revmap[bphf->lookup(unikmer_hash)] = unikmer_hash;
    }
    //std::cerr << "Finished building MPHF." << std::endl;
}

UKHS::Map::~Map() = default;

bool UKHS::Map::query(Unikmer& unikmer) {
    unikmer.partition = ULLONG_MAX;
    unikmer.partition = bphf->lookup(unikmer.hash);
    if (!unikmer.is_valid()) {
        return false;
    }
    if (ukhs_revmap[unikmer.partition] == unikmer.hash) {
        return true;
    }
    return false;
}


}
}
