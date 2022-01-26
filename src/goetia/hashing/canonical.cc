/**
 * (c) Camille Scott, 2019
 * File   : canonical.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 08.01.2020
 */

#include "goetia/hashing/canonical.hh"



namespace goetia {

    template class goetia::Hash<uint64_t>;
    template class goetia::Canonical<goetia::Hash<uint64_t>>;

    template class goetia::Kmer<goetia::Hash<uint64_t>>;
    template class goetia::Kmer<goetia::Canonical<uint64_t>>;

    template class goetia::Kmer<goetia::UnikmerWmer>;
    template class goetia::Kmer<goetia::CanonicalUnikmerWmer>;

    template class goetia::Wmer<goetia::Hash<uint64_t>, goetia::Unikmer>;
    template class goetia::Wmer<goetia::Canonical<uint64_t>, goetia::CanonicalUnikmer>;

    template class goetia::Shift<goetia::Hash<uint64_t>, goetia::DIR_LEFT>;
    template class goetia::Shift<goetia::Hash<uint64_t>, goetia::DIR_RIGHT>;
    template class goetia::Shift<goetia::Canonical<uint64_t>, goetia::DIR_LEFT>;
    template class goetia::Shift<goetia::Canonical<uint64_t>, goetia::DIR_RIGHT>;

}

