/**
 * (c) Camille Scott, 2019
 * File   : saturating_compactor.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.10.2019
 */

#include "boink/cdbg/saturating_compactor.hh"

#include "boink/dbg.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/storage/storage_types.hh"
#include "boink/signatures/sourmash_signature.hh"


template class boink::cdbg::SaturatingCompactor<boink::dBG<boink::storage::BitStorage,
                                                           boink::hashing::RollingHashShifter>,
                                                boink::signatures::SourmashSignature::Signature>;
template class boink::cdbg::SaturatingCompactor<boink::dBG<boink::storage::ByteStorage,
                                                           boink::hashing::RollingHashShifter>,
                                                boink::signatures::SourmashSignature::Signature>;
template class boink::cdbg::SaturatingCompactor<boink::dBG<boink::storage::NibbleStorage,
                                                           boink::hashing::RollingHashShifter>,
                                                boink::signatures::SourmashSignature::Signature>;
template class boink::cdbg::SaturatingCompactor<boink::dBG<boink::storage::QFStorage,
                                                           boink::hashing::RollingHashShifter>,
                                                boink::signatures::SourmashSignature::Signature>;
template class boink::cdbg::SaturatingCompactor<boink::dBG<boink::storage::SparseppSetStorage,
                                                           boink::hashing::RollingHashShifter>,
                                                boink::signatures::SourmashSignature::Signature>;
