/* sparseppstorage.hh -- storage classes for the boink dbg
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_SPARSEPPSETSTORAGE_HH
#define BOINK_SPARSEPPSETSTORAGE_HH

#include "boink/hashing/hashing_types.hh"
#include "boink/storage/storage.hh"
#include "sparsepp/spp.h"

#include <vector>

namespace boink {
namespace storage {


class SparseppSetStorage : public Storage {

protected:

    spp::sparse_hash_set<hashing::hash_t> _store;

public:
    
    SparseppSetStorage(const std::vector<uint64_t>& sizes)
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

    void save(std::string, uint16_t ) {
     
    }

    void load(std::string, uint16_t &) {

    }

    count_t test_and_set_bits(hashing::hash_t h) {
        count_t count = get_count(h);
        add(h);
        return count;
    }

    bool add(hashing::hash_t h) {
        count_t count = get_count(h);
        _store.insert(h);
        return count == 0;
    }

    const count_t get_count(hashing::hash_t h) const {
        return _store.count(h);
    }

    byte_t ** get_raw_tables() {
        return nullptr;
    }

};

}
}
#endif
