/* partitioned_storage.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "goetia/storage/partitioned_storage.hh"
#include "goetia/storage/storage_types.hh"


template class goetia::storage::PartitionedStorage<goetia::storage::BitStorage>;
template class goetia::storage::PartitionedStorage<goetia::storage::ByteStorage>;
template class goetia::storage::PartitionedStorage<goetia::storage::NibbleStorage>;
template class goetia::storage::PartitionedStorage<goetia::storage::QFStorage>;
template class goetia::storage::PartitionedStorage<goetia::storage::SparseppSetStorage>;
template class goetia::storage::PartitionedStorage<goetia::storage::PHMapStorage>;
template class goetia::storage::PartitionedStorage<goetia::storage::BTreeStorage>;
