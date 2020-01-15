/* compactor.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "boink/cdbg/compactor.hh"

#include "boink/dbg.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/storage/storage_types.hh"

namespace boink {

template class cdbg::StreamingCompactor<dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>>;
//template class cdbg::StreamingCompactor<dBG<storage::SparseppSetStorage, hashing::CanRollingShifter>>;

template class cdbg::StreamingCompactor<dBG<storage::BitStorage, hashing::FwdRollingShifter>>;
//template class cdbg::StreamingCompactor<dBG<storage::BitStorage, hashing::CanRollingShifter>>;

template class cdbg::StreamingCompactor<dBG<storage::ByteStorage, hashing::FwdRollingShifter>>;
//template class cdbg::StreamingCompactor<dBG<storage::ByteStorage, hashing::CanRollingShifter>>;

template class cdbg::StreamingCompactor<dBG<storage::NibbleStorage, hashing::FwdRollingShifter>>;
//template class cdbg::StreamingCompactor<dBG<storage::NibbleStorage, hashing::CanRollingShifter>>;

template class cdbg::StreamingCompactor<dBG<storage::QFStorage, hashing::FwdRollingShifter>>;
//template class cdbg::StreamingCompactor<dBG<storage::QFStorage, hashing::CanRollingShifter>>;
}
