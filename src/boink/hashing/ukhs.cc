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

UKHS::UKHS(uint16_t K,
           std::vector<std::string>& ukhs)
    : KmerClient  (K),
      ukhs_revmap (ukhs.size()),
      ukhs_hasher (K)
{
    if (ukhs.front().size() != K) {
        throw BoinkException("K does not match k-mer size from provided UKHS");
    }

    //std::cerr << "Building MPHF with "
    //          << ukhs.size() << " k-mers."
    //          << std::endl;

    for (auto unikmer : ukhs) {
        ukhs_hashes.push_back(ukhs_hasher.hash(unikmer));
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

UKHS::~UKHS() = default;

bool UKHS::query(Unikmer& unikmer) {
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
