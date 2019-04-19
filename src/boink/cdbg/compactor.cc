/* boink.hh
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
#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/sparseppstorage.hh"


template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::BitStorage,
                                               boink::hashing::RollingHashShifter>>;
template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::ByteStorage,
                                               boink::hashing::RollingHashShifter>>;
template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::NibbleStorage,
                                               boink::hashing::RollingHashShifter>>;
template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::QFStorage,
                                               boink::hashing::RollingHashShifter>>;
template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::SparseppSetStorage,
                                               boink::hashing::RollingHashShifter>>;
