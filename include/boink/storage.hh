/* storage.hh -- storage classes for the boink dbg
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_STORAGE_HH
#define BOINK_STORAGE_HH

#include "hashing.hh"
#include "oxli/storage.hh"
#include "sparsepp/sparsepp/spp.h"

#include <vector>

namespace boink {

typedef oxli::BoundedCounterType count_t;
typedef std::pair<uint8_t, uint8_t> full_count_t;


class SparseppSetStorage : public oxli::Storage {

protected:

    spp::sparse_hash_set<hash_t> _store;

public:
    
    SparseppSetStorage(std::vector<uint64_t>& sizes)
    {
    }

    void reset() {
        _store.clear();
    }

    std::vector<uint64_t> get_tablesizes() const {
        return std::vector<uint64_t>({_store.max_size()});
    }

    const uint64_t n_unique_kmers() const {
        return _store.size();
    }

    const uint64_t n_tables() const {
        return 1;
    }

    const uint64_t n_occupied() const {
        return _store.bucket_count();
    }

    void save(std::string, oxli::WordLength) {
     
    }

    void load(std::string, oxli::WordLength&) {

    }

    count_t test_and_set_bits(oxli::HashIntoType h) {
        count_t count = get_count(h);
        add(h);
        return count;
    }

    bool add(oxli::HashIntoType h) {
        count_t count = get_count(h);
        _store.insert(h);
        return count == 0;
    }

    const count_t get_count(oxli::HashIntoType h) const {
        return _store.count(h);
    }

    oxli::Byte ** get_raw_tables() {
        return nullptr;
    }

};

}
#endif
