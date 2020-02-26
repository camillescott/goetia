/**
 * (c) Camille Scott, 2019
 * File   : cdbg.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 03.09.2019
 */

#include "goetia/cdbg/cdbg.hh"

#include "goetia/dbg.hh"
#include "goetia/hashing/rollinghashshifter.hh"
#include "goetia/storage/storage_types.hh"

namespace goetia {

template class cdbg::cDBG<dBG<storage::BitStorage, hashing::FwdLemireShifter>>;
template class cdbg::cDBG<dBG<storage::BitStorage, hashing::CanLemireShifter>>;

template class cdbg::cDBG<dBG<storage::SparseppSetStorage, hashing::FwdLemireShifter>>;
template class cdbg::cDBG<dBG<storage::SparseppSetStorage, hashing::CanLemireShifter>>;

template class cdbg::cDBG<dBG<storage::ByteStorage, hashing::FwdLemireShifter>>;
template class cdbg::cDBG<dBG<storage::ByteStorage, hashing::CanLemireShifter>>;

template class cdbg::cDBG<dBG<storage::NibbleStorage, hashing::FwdLemireShifter>>;
template class cdbg::cDBG<dBG<storage::NibbleStorage, hashing::CanLemireShifter>>;

template class cdbg::cDBG<dBG<storage::QFStorage, hashing::FwdLemireShifter>>;
template class cdbg::cDBG<dBG<storage::QFStorage, hashing::CanLemireShifter>>;
}
