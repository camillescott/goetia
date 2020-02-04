/**
 * (c) Camille Scott, 2019
 * File   : cdbg.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 03.09.2019
 */

#include "boink/cdbg/cdbg.hh"

#include "boink/dbg.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/storage/storage_types.hh"

namespace boink {

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
