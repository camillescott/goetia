/**
 * (c) Camille Scott, 2019
 * File   : canonical.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 08.01.2020
 */

#include "boink/hashing/canonical.hh"



namespace boink::hashing {

    template class boink::hashing::Hash<uint64_t>;
    template class boink::hashing::Canonical<boink::hashing::Hash<uint64_t>>;

    template class boink::hashing::Kmer<boink::hashing::Hash<uint64_t>>;
    template class boink::hashing::Kmer<boink::hashing::Canonical<uint64_t>>;

    template class boink::hashing::Kmer<boink::hashing::UnikmerWmer>;
    template class boink::hashing::Kmer<boink::hashing::CanonicalUnikmerWmer>;

    template class boink::hashing::Wmer<boink::hashing::Hash<uint64_t>, boink::hashing::Unikmer>;
    template class boink::hashing::Wmer<boink::hashing::Canonical<uint64_t>, boink::hashing::CanonicalUnikmer>;

    template class boink::hashing::Shift<boink::hashing::Hash<uint64_t>, boink::hashing::DIR_LEFT>;
    template class boink::hashing::Shift<boink::hashing::Hash<uint64_t>, boink::hashing::DIR_RIGHT>;
    template class boink::hashing::Shift<boink::hashing::Canonical<uint64_t>, boink::hashing::DIR_LEFT>;
    template class boink::hashing::Shift<boink::hashing::Canonical<uint64_t>, boink::hashing::DIR_RIGHT>;

}

