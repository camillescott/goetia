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

    template class dBG<storage::BitStorage, hashing::FwdLemireShifter>;
    template class Masked<storage::BitStorage, hashing::FwdLemireShifter, std::set<typename hashing::FwdLemireShifter::hash_type>>;
    template class dBG<storage::BitStorage, hashing::CanLemireShifter>;
    template class dBG<storage::BitStorage, hashing::FwdUnikmerShifter>;
    template class dBG<storage::BitStorage, hashing::CanUnikmerShifter>;

    template class dBG<storage::SparseppSetStorage, hashing::FwdLemireShifter>;
    template class dBG<storage::SparseppSetStorage, hashing::CanLemireShifter>;
    template class dBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>;
    template class dBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>;

    template class dBG<storage::ByteStorage, hashing::FwdLemireShifter>;
    template class dBG<storage::ByteStorage, hashing::CanLemireShifter>;
    template class dBG<storage::ByteStorage, hashing::FwdUnikmerShifter>;
    template class dBG<storage::ByteStorage, hashing::CanUnikmerShifter>;

    template class dBG<storage::NibbleStorage, hashing::FwdLemireShifter>;
    template class dBG<storage::NibbleStorage, hashing::CanLemireShifter>;
    template class dBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>;
    template class dBG<storage::NibbleStorage, hashing::CanUnikmerShifter>;

    template class dBG<storage::QFStorage, hashing::FwdLemireShifter>;
    template class dBG<storage::QFStorage, hashing::CanLemireShifter>;
    template class dBG<storage::QFStorage, hashing::FwdUnikmerShifter>;
    template class dBG<storage::QFStorage, hashing::CanUnikmerShifter>;

    template class goetia::dBG<goetia::storage::PHMapStorage, goetia::hashing::FwdLemireShifter>;
    template class goetia::dBG<goetia::storage::PHMapStorage, goetia::hashing::CanLemireShifter>;
    template class goetia::dBG<goetia::storage::PHMapStorage, goetia::hashing::FwdUnikmerShifter>;
    template class goetia::dBG<goetia::storage::PHMapStorage, goetia::hashing::CanUnikmerShifter>;
   
    template class goetia::dBG<goetia::storage::BTreeStorage, goetia::hashing::FwdLemireShifter>;
    template class goetia::dBG<goetia::storage::BTreeStorage, goetia::hashing::CanLemireShifter>;
    template class goetia::dBG<goetia::storage::BTreeStorage, goetia::hashing::FwdUnikmerShifter>;
    template class goetia::dBG<goetia::storage::BTreeStorage, goetia::hashing::CanUnikmerShifter>;


    template class hashing::KmerIterator<dBG<storage::BitStorage, hashing::FwdLemireShifter>>;
    template class hashing::KmerIterator<dBG<storage::BitStorage, hashing::CanLemireShifter>>;
    template class hashing::KmerIterator<dBG<storage::BitStorage, hashing::FwdUnikmerShifter>>;
    template class hashing::KmerIterator<dBG<storage::BitStorage, hashing::CanUnikmerShifter>>;

    template class hashing::KmerIterator<dBG<storage::SparseppSetStorage, hashing::FwdLemireShifter>>;
    template class hashing::KmerIterator<dBG<storage::SparseppSetStorage, hashing::CanLemireShifter>>;
    template class hashing::KmerIterator<dBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>>;
    template class hashing::KmerIterator<dBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>>;

    template class hashing::KmerIterator<dBG<storage::ByteStorage, hashing::FwdLemireShifter>>;
    template class hashing::KmerIterator<dBG<storage::ByteStorage, hashing::CanLemireShifter>>;
    template class hashing::KmerIterator<dBG<storage::ByteStorage, hashing::FwdUnikmerShifter>>;
    template class hashing::KmerIterator<dBG<storage::ByteStorage, hashing::CanUnikmerShifter>>;

    template class hashing::KmerIterator<dBG<storage::NibbleStorage, hashing::FwdLemireShifter>>;
    template class hashing::KmerIterator<dBG<storage::NibbleStorage, hashing::CanLemireShifter>>;
    template class hashing::KmerIterator<dBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>>;
    template class hashing::KmerIterator<dBG<storage::NibbleStorage, hashing::CanUnikmerShifter>>;

    template class hashing::KmerIterator<dBG<storage::QFStorage, hashing::FwdLemireShifter>>;
    template class hashing::KmerIterator<dBG<storage::QFStorage, hashing::CanLemireShifter>>;
    template class hashing::KmerIterator<dBG<storage::QFStorage, hashing::FwdUnikmerShifter>>;
    template class hashing::KmerIterator<dBG<storage::QFStorage, hashing::CanUnikmerShifter>>;

    template class goetia::hashing::KmerIterator<goetia::dBG<goetia::storage::PHMapStorage, goetia::hashing::FwdLemireShifter>>;
    template class goetia::hashing::KmerIterator<goetia::dBG<goetia::storage::PHMapStorage, goetia::hashing::CanLemireShifter>>;
    template class goetia::hashing::KmerIterator<goetia::dBG<goetia::storage::PHMapStorage, goetia::hashing::FwdUnikmerShifter>>;
    template class goetia::hashing::KmerIterator<goetia::dBG<goetia::storage::PHMapStorage, goetia::hashing::CanUnikmerShifter>>;

    template class goetia::hashing::KmerIterator<goetia::dBG<goetia::storage::BTreeStorage, goetia::hashing::FwdLemireShifter>>;
    template class goetia::hashing::KmerIterator<goetia::dBG<goetia::storage::BTreeStorage, goetia::hashing::CanLemireShifter>>;
    template class goetia::hashing::KmerIterator<goetia::dBG<goetia::storage::BTreeStorage, goetia::hashing::FwdUnikmerShifter>>;
    template class goetia::hashing::KmerIterator<goetia::dBG<goetia::storage::BTreeStorage, goetia::hashing::CanUnikmerShifter>>;


}
