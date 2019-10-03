/* boink/cdbg/ucompactor.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "boink/cdbg/ucompactor.hh"

#include "boink/dbg.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/storage/storage_types.hh"

template class boink::cdbg::UStreamingCompactor<boink::storage::BitStorage<>>;
template class boink::cdbg::UStreamingCompactor<boink::storage::ByteStorage>;
template class boink::cdbg::UStreamingCompactor<boink::storage::NibbleStorage<>>;
template class boink::cdbg::UStreamingCompactor<boink::storage::QFStorage<>>;
template class boink::cdbg::UStreamingCompactor<boink::storage::SparseppSetStorage<>>;
