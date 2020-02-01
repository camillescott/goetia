/**
 * (c) Camille Scott, 2019
 * File   : compactor.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.01.2020
 */

#include "boink/cdbg/compactor.hh"
#include "boink/storage/storage_types.hh"
#include "boink/hashing/rollinghashshifter.hh"


template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::SparseppSetStorage, boink::hashing::FwdLemireShifter>>;
template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::BitStorage, boink::hashing::FwdLemireShifter>>;
template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::ByteStorage, boink::hashing::FwdLemireShifter>>;
template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::NibbleStorage, boink::hashing::FwdLemireShifter>>;
template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::QFStorage, boink::hashing::FwdLemireShifter>>;

