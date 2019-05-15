/* boink.hh
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "boink/storage/sparseppstorage.hh"
#include "boink/storage/sparsepp/spp.h"

#include <cstdint>

using namespace boink;
using namespace boink::hashing;
using namespace boink::storage;

const bool
SparseppSetStorage::insert(hash_t h) {
    auto result = _store->insert(h);
    // the second in the returned pair reports that the insert
    // took place ie the hash was new
    return result.second;
}


const count_t
SparseppSetStorage::insert_and_query(hashing::hash_t h) {
    insert(h);
    return 1; // its a presence filter so always 1 after insert
}


const count_t
SparseppSetStorage::query(hashing::hash_t h) const {
    return _store->count(h);
}
