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


template class goetia::PartitionedStorage<goetia::BitStorage>;
template class goetia::PartitionedStorage<goetia::ByteStorage>;
template class goetia::PartitionedStorage<goetia::NibbleStorage>;
template class goetia::PartitionedStorage<goetia::QFStorage>;
template class goetia::PartitionedStorage<goetia::SparseppSetStorage>;
template class goetia::PartitionedStorage<goetia::PHMapStorage>;
template class goetia::PartitionedStorage<goetia::BTreeStorage>;
