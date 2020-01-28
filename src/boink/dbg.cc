/**
 * (c) Camille Scott, 2019
 * File   : dbg.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */

#include "boink/dbg.hh"

#include "boink/storage/storage_types.hh"
#include "boink/hashing/shifter_types.hh"


namespace boink {

    template class dBG<storage::BitStorage, hashing::FwdRollingShifter>;
    template class dBG<storage::BitStorage, hashing::CanRollingShifter>;
    template class dBG<storage::BitStorage, hashing::FwdUnikmerShifter>;
    template class dBG<storage::BitStorage, hashing::CanUnikmerShifter>;

    template class dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>;
    template class dBG<storage::SparseppSetStorage, hashing::CanRollingShifter>;
    template class dBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>;
    template class dBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>;

    template class dBG<storage::ByteStorage, hashing::FwdRollingShifter>;
    template class dBG<storage::ByteStorage, hashing::CanRollingShifter>;
    template class dBG<storage::ByteStorage, hashing::FwdUnikmerShifter>;
    template class dBG<storage::ByteStorage, hashing::CanUnikmerShifter>;

    template class dBG<storage::NibbleStorage, hashing::FwdRollingShifter>;
    template class dBG<storage::NibbleStorage, hashing::CanRollingShifter>;
    template class dBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>;
    template class dBG<storage::NibbleStorage, hashing::CanUnikmerShifter>;

    template class dBG<storage::QFStorage, hashing::FwdRollingShifter>;
    template class dBG<storage::QFStorage, hashing::CanRollingShifter>;
    template class dBG<storage::QFStorage, hashing::FwdUnikmerShifter>;
    template class dBG<storage::QFStorage, hashing::CanUnikmerShifter>;


    template class hashing::KmerIterator<dBG<storage::BitStorage, hashing::FwdRollingShifter>>;
    template class hashing::KmerIterator<dBG<storage::BitStorage, hashing::CanRollingShifter>>;
    template class hashing::KmerIterator<dBG<storage::BitStorage, hashing::FwdUnikmerShifter>>;
    template class hashing::KmerIterator<dBG<storage::BitStorage, hashing::CanUnikmerShifter>>;

    template class hashing::KmerIterator<dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>>;
    template class hashing::KmerIterator<dBG<storage::SparseppSetStorage, hashing::CanRollingShifter>>;
    template class hashing::KmerIterator<dBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>>;
    template class hashing::KmerIterator<dBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>>;

    template class hashing::KmerIterator<dBG<storage::ByteStorage, hashing::FwdRollingShifter>>;
    template class hashing::KmerIterator<dBG<storage::ByteStorage, hashing::CanRollingShifter>>;
    template class hashing::KmerIterator<dBG<storage::ByteStorage, hashing::FwdUnikmerShifter>>;
    template class hashing::KmerIterator<dBG<storage::ByteStorage, hashing::CanUnikmerShifter>>;

    template class hashing::KmerIterator<dBG<storage::NibbleStorage, hashing::FwdRollingShifter>>;
    template class hashing::KmerIterator<dBG<storage::NibbleStorage, hashing::CanRollingShifter>>;
    template class hashing::KmerIterator<dBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>>;
    template class hashing::KmerIterator<dBG<storage::NibbleStorage, hashing::CanUnikmerShifter>>;

    template class hashing::KmerIterator<dBG<storage::QFStorage, hashing::FwdRollingShifter>>;
    template class hashing::KmerIterator<dBG<storage::QFStorage, hashing::CanRollingShifter>>;
    template class hashing::KmerIterator<dBG<storage::QFStorage, hashing::FwdUnikmerShifter>>;
    template class hashing::KmerIterator<dBG<storage::QFStorage, hashing::CanUnikmerShifter>>;
}
