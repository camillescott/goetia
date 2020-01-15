/* boink/cdbg/udbg.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "boink/cdbg/udbg.hh"

#include "boink/dbg.hh"
#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/sparseppstorage.hh"

#include "boink/hashing/rollinghashshifter.hh"

#include "boink/boink.hh"

namespace boink {

template class cdbg::uDBG<storage::BitStorage>;
//template class cdbg::uDBG<storage::ByteStorage>;
//template class cdbg::uDBG<storage::NibbleStorage>;
//template class cdbg::uDBG<storage::QFStorage>;
//template class cdbg::uDBG<storage::SparseppSetStorage>;

}




