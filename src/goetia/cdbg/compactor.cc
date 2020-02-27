/**
 * (c) Camille Scott, 2019
 * File   : compactor.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.01.2020
 */

#include "goetia/cdbg/compactor.hh"
#include "goetia/storage/storage_types.hh"
#include "goetia/hashing/rollinghashshifter.hh"


template class goetia::cdbg::StreamingCompactor<goetia::dBG<goetia::storage::SparseppSetStorage, goetia::hashing::FwdLemireShifter>>;
// template class goetia::cdbg::StreamingCompactor<goetia::dBG<goetia::storage::BitStorage, goetia::hashing::FwdLemireShifter>>;
// template class goetia::cdbg::StreamingCompactor<goetia::dBG<goetia::storage::ByteStorage, goetia::hashing::FwdLemireShifter>>;
// template class goetia::cdbg::StreamingCompactor<goetia::dBG<goetia::storage::NibbleStorage, goetia::hashing::FwdLemireShifter>>;
// template class goetia::cdbg::StreamingCompactor<goetia::dBG<goetia::storage::QFStorage, goetia::hashing::FwdLemireShifter>>;

