/* ukhs_signature.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "boink/signatures/ukhs_signature.hh"

#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/partitioned_storage.hh"
#include "boink/storage/sparseppstorage.hh"


using namespace boink;
using namespace boink::signatures;


template class boink::signatures::UnikmerSignature<boink::storage::BitStorage>;
template class boink::signatures::UnikmerSignature<boink::storage::ByteStorage>;
template class boink::signatures::UnikmerSignature<boink::storage::NibbleStorage>;
template class boink::signatures::UnikmerSignature<boink::storage::QFStorage>;
template class boink::signatures::UnikmerSignature<boink::storage::SparseppSetStorage>;
