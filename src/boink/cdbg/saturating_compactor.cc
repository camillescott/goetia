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


namespace boink {

template class cdbg::SaturatingCompactor<dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>,
                                         signatures::SourmashSignature::Signature>;

}

