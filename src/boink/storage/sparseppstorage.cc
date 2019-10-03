/**
 * (c) Camille Scott, 2019
 * File   : sparseppstorage.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */
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

namespace boink {
namespace storage {

template<class ValueType>
const bool
SparseppSetStorage<ValueType>::insert(value_type h) {
    auto result = _store->insert(h);
    // the second in the returned pair reports that the insert
    // took place ie the hash was new
    return result.second;
}


template<class ValueType>
const count_t
SparseppSetStorage<ValueType>::insert_and_query(value_type h) {
    insert(h);
    return 1; // its a presence filter so always 1 after insert
}


template<class ValueType>
const count_t
SparseppSetStorage<ValueType>::query(value_type h) const {
    return _store->count(h);
}


template<class ValueType>
std::shared_ptr<SparseppSetStorage<ValueType>>
SparseppSetStorage<ValueType>::build() {
    return std::make_shared<SparseppSetStorage>();
}


template<class ValueType>
std::shared_ptr<SparseppSetStorage<ValueType>>
SparseppSetStorage<ValueType>::clone() const {
    return std::make_shared<SparseppSetStorage>();
}

}
}

template class boink::storage::SparseppSetStorage<uint64_t>;
