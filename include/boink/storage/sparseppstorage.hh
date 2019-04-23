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

#include "boink/storage/storage.hh"
#include "boink/storage/sparsepp/spp.h"

#include <cstdint>
#include <memory>
#include <vector>

namespace boink {
namespace storage {


class SparseppSetStorage : public Storage {

public:

    typedef spp::sparse_hash_set<hashing::hash_t> store_type;

protected:

    std::unique_ptr<store_type> _store;

public:
    
    SparseppSetStorage()
    {
        _store = std::make_unique<store_type>();
    }

    std::unique_ptr<SparseppSetStorage> clone() const {
        return std::make_unique<SparseppSetStorage>();
    }

    void reset() {
        _store->clear();
    }

    const uint64_t get_maxsize() const {
        return _store->max_size();
    }

    const uint64_t n_unique_kmers() const {
        return _store->size();
    }

    const uint64_t n_buckets() const {
        return _store->bucket_count();
    }

    const uint64_t n_occupied() const {
        return n_buckets();
    }

    void save(std::string, uint16_t ) {
     
    }

    void load(std::string, uint16_t &) {

    }

    const bool insert(hashing::hash_t h) {
        auto result = _store->insert(h);
        // the second in the returned pair reports that the insert
        // took place ie the hash was new
        return result.second;
    }

    const count_t insert_and_query(hashing::hash_t h) {
        insert(h);
        return 1; // its a presence filter so always 1 after insert
    }

    const count_t query(hashing::hash_t h) const {
        return _store->count(h);
    }


    byte_t ** get_raw_tables() {
        return nullptr;
    }

};


template<> 
struct is_probabilistic<SparseppSetStorage> { 
      static const bool value = false;
};


}
}
#endif
