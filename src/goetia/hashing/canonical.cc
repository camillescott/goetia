/**
 * (c) Camille Scott, 2019
 * File   : canonical.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 08.01.2020
 */

#include "goetia/hashing/canonical.hh"



namespace goetia::hashing {

    template class goetia::hashing::Hash<uint64_t>;
    template class goetia::hashing::Canonical<goetia::hashing::Hash<uint64_t>>;

    template class goetia::hashing::Kmer<goetia::hashing::Hash<uint64_t>>;
    template class goetia::hashing::Kmer<goetia::hashing::Canonical<uint64_t>>;

    template class goetia::hashing::Kmer<goetia::hashing::UnikmerWmer>;
    template class goetia::hashing::Kmer<goetia::hashing::CanonicalUnikmerWmer>;

    template class goetia::hashing::Wmer<goetia::hashing::Hash<uint64_t>, goetia::hashing::Unikmer>;
    template class goetia::hashing::Wmer<goetia::hashing::Canonical<uint64_t>, goetia::hashing::CanonicalUnikmer>;

    template class goetia::hashing::Shift<goetia::hashing::Hash<uint64_t>, goetia::hashing::DIR_LEFT>;
    template class goetia::hashing::Shift<goetia::hashing::Hash<uint64_t>, goetia::hashing::DIR_RIGHT>;
    template class goetia::hashing::Shift<goetia::hashing::Canonical<uint64_t>, goetia::hashing::DIR_LEFT>;
    template class goetia::hashing::Shift<goetia::hashing::Canonical<uint64_t>, goetia::hashing::DIR_RIGHT>;

}

