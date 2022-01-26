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


template class goetia::StreamingCompactor<goetia::dBG<goetia::SparseppSetStorage, goetia::FwdLemireShifter>>;
template class goetia::StreamingCompactor<goetia::dBG<goetia::PHMapStorage, goetia::FwdLemireShifter>>;
// template class goetia::StreamingCompactor<goetia::dBG<goetia::BitStorage, goetia::FwdLemireShifter>>;
// template class goetia::StreamingCompactor<goetia::dBG<goetia::ByteStorage, goetia::FwdLemireShifter>>;
// template class goetia::StreamingCompactor<goetia::dBG<goetia::NibbleStorage, goetia::FwdLemireShifter>>;
// template class goetia::StreamingCompactor<goetia::dBG<goetia::QFStorage, goetia::FwdLemireShifter>>;

