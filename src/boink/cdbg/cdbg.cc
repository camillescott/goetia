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

template class cdbg::cDBG<dBG<storage::BitStorage, hashing::FwdRollingShifter>>;
template class cdbg::cDBG<dBG<storage::BitStorage, hashing::CanRollingShifter>>;

template class cdbg::cDBG<dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>>;
template class cdbg::cDBG<dBG<storage::SparseppSetStorage, hashing::CanRollingShifter>>;

template class cdbg::cDBG<dBG<storage::ByteStorage, hashing::FwdRollingShifter>>;
template class cdbg::cDBG<dBG<storage::ByteStorage, hashing::CanRollingShifter>>;

template class cdbg::cDBG<dBG<storage::NibbleStorage, hashing::FwdRollingShifter>>;
template class cdbg::cDBG<dBG<storage::NibbleStorage, hashing::CanRollingShifter>>;

template class cdbg::cDBG<dBG<storage::QFStorage, hashing::FwdRollingShifter>>;
template class cdbg::cDBG<dBG<storage::QFStorage, hashing::CanRollingShifter>>;
}
