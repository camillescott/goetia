/**
 * (c) Camille Scott, 2019
 * File   : dbg.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */

#include "goetia/dbg.hh"

#include "goetia/storage/storage_types.hh"
#include "goetia/hashing/shifter_types.hh"


namespace goetia {

    template class dBG<BitStorage, FwdLemireShifter>;
    template class Masked<BitStorage, FwdLemireShifter, std::set<typename FwdLemireShifter::hash_type>>;
    template class dBG<BitStorage, CanLemireShifter>;
    template class dBG<BitStorage, FwdUnikmerShifter>;
    template class dBG<BitStorage, CanUnikmerShifter>;

    template class dBG<SparseppSetStorage, FwdLemireShifter>;
    template class dBG<SparseppSetStorage, CanLemireShifter>;
    template class dBG<SparseppSetStorage, FwdUnikmerShifter>;
    template class dBG<SparseppSetStorage, CanUnikmerShifter>;

    template class dBG<ByteStorage, FwdLemireShifter>;
    template class dBG<ByteStorage, CanLemireShifter>;
    template class dBG<ByteStorage, FwdUnikmerShifter>;
    template class dBG<ByteStorage, CanUnikmerShifter>;

    template class dBG<NibbleStorage, FwdLemireShifter>;
    template class dBG<NibbleStorage, CanLemireShifter>;
    template class dBG<NibbleStorage, FwdUnikmerShifter>;
    template class dBG<NibbleStorage, CanUnikmerShifter>;

    template class dBG<QFStorage, FwdLemireShifter>;
    template class dBG<QFStorage, CanLemireShifter>;
    template class dBG<QFStorage, FwdUnikmerShifter>;
    template class dBG<QFStorage, CanUnikmerShifter>;

    template class goetia::dBG<goetia::PHMapStorage, goetia::FwdLemireShifter>;
    template class goetia::dBG<goetia::PHMapStorage, goetia::CanLemireShifter>;
    template class goetia::dBG<goetia::PHMapStorage, goetia::FwdUnikmerShifter>;
    template class goetia::dBG<goetia::PHMapStorage, goetia::CanUnikmerShifter>;
   
    template class goetia::dBG<goetia::BTreeStorage, goetia::FwdLemireShifter>;
    template class goetia::dBG<goetia::BTreeStorage, goetia::CanLemireShifter>;
    template class goetia::dBG<goetia::BTreeStorage, goetia::FwdUnikmerShifter>;
    template class goetia::dBG<goetia::BTreeStorage, goetia::CanUnikmerShifter>;


    template class KmerIterator<dBG<BitStorage, FwdLemireShifter>>;
    template class KmerIterator<dBG<BitStorage, CanLemireShifter>>;
    template class KmerIterator<dBG<BitStorage, FwdUnikmerShifter>>;
    template class KmerIterator<dBG<BitStorage, CanUnikmerShifter>>;

    template class KmerIterator<dBG<SparseppSetStorage, FwdLemireShifter>>;
    template class KmerIterator<dBG<SparseppSetStorage, CanLemireShifter>>;
    template class KmerIterator<dBG<SparseppSetStorage, FwdUnikmerShifter>>;
    template class KmerIterator<dBG<SparseppSetStorage, CanUnikmerShifter>>;

    template class KmerIterator<dBG<ByteStorage, FwdLemireShifter>>;
    template class KmerIterator<dBG<ByteStorage, CanLemireShifter>>;
    template class KmerIterator<dBG<ByteStorage, FwdUnikmerShifter>>;
    template class KmerIterator<dBG<ByteStorage, CanUnikmerShifter>>;

    template class KmerIterator<dBG<NibbleStorage, FwdLemireShifter>>;
    template class KmerIterator<dBG<NibbleStorage, CanLemireShifter>>;
    template class KmerIterator<dBG<NibbleStorage, FwdUnikmerShifter>>;
    template class KmerIterator<dBG<NibbleStorage, CanUnikmerShifter>>;

    template class KmerIterator<dBG<QFStorage, FwdLemireShifter>>;
    template class KmerIterator<dBG<QFStorage, CanLemireShifter>>;
    template class KmerIterator<dBG<QFStorage, FwdUnikmerShifter>>;
    template class KmerIterator<dBG<QFStorage, CanUnikmerShifter>>;

    template class goetia::KmerIterator<goetia::dBG<goetia::PHMapStorage, goetia::FwdLemireShifter>>;
    template class goetia::KmerIterator<goetia::dBG<goetia::PHMapStorage, goetia::CanLemireShifter>>;
    template class goetia::KmerIterator<goetia::dBG<goetia::PHMapStorage, goetia::FwdUnikmerShifter>>;
    template class goetia::KmerIterator<goetia::dBG<goetia::PHMapStorage, goetia::CanUnikmerShifter>>;

    template class goetia::KmerIterator<goetia::dBG<goetia::BTreeStorage, goetia::FwdLemireShifter>>;
    template class goetia::KmerIterator<goetia::dBG<goetia::BTreeStorage, goetia::CanLemireShifter>>;
    template class goetia::KmerIterator<goetia::dBG<goetia::BTreeStorage, goetia::FwdUnikmerShifter>>;
    template class goetia::KmerIterator<goetia::dBG<goetia::BTreeStorage, goetia::CanUnikmerShifter>>;


}
