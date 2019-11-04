/**
 * (c) Camille Scott, 2019
 * File   : sparseppstorage.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */

#include "boink/storage/sparseppstorage.hh"
#include "boink/storage/sparsepp/spp.h"

#include <cstdint>

namespace boink {
namespace storage {

const bool
SparseppSetStorage::insert(value_type h) {
    auto result = _store->insert(h);
    // the second in the returned pair reports that the insert
    // took place ie the hash was new
    return result.second;
}


const count_t
SparseppSetStorage::insert_and_query(value_type h) {
    insert(h);
    return 1; // its a presence filter so always 1 after insert
}


const count_t
SparseppSetStorage::query(value_type h) const {
    return _store->count(h);
}


std::shared_ptr<SparseppSetStorage>
SparseppSetStorage::build() {
    return std::make_shared<SparseppSetStorage>();
}


std::shared_ptr<SparseppSetStorage>
SparseppSetStorage::clone() const {
    return std::make_shared<SparseppSetStorage>();
}

}
}
