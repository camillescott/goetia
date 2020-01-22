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

namespace boink {

    template class cdbg::StreamingCompactor<dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>>;
    template class cdbg::StreamingCompactor<dBG<storage::BitStorage, hashing::FwdRollingShifter>>;
    template class cdbg::StreamingCompactor<dBG<storage::ByteStorage, hashing::FwdRollingShifter>>;
    template class cdbg::StreamingCompactor<dBG<storage::NibbleStorage, hashing::FwdRollingShifter>>;
    template class cdbg::StreamingCompactor<dBG<storage::QFStorage, hashing::FwdRollingShifter>>;

}
