/**
 * (c) Camille Scott, 2019
 * File   : sparseppstorage.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */


#ifndef GOETIA_SPARSEPPSETSTORAGE_HH
#define GOETIA_SPARSEPPSETSTORAGE_HH

#include "goetia/meta.hh"
#include "goetia/storage/storage.hh"
#include "goetia/storage/sparsepp/spp.h"

#include <cstdint>
#include <fstream>
#include <memory>
#include <tuple>
#include <vector>

namespace goetia::storage {

class SparseppSetStorage;


template<>
struct StorageTraits<SparseppSetStorage> {
    static constexpr bool is_probabilistic = false;
    static constexpr bool is_counting      = false;

    typedef std::tuple<bool> params_type;
    static constexpr params_type default_params = std::make_tuple(0);
};


class SparseppSetStorage : public Storage<uint64_t>,
                           public Tagged<SparseppSetStorage> {

public:
    
    using Storage<uint64_t>::value_type;
    using Traits = StorageTraits<SparseppSetStorage>;
    typedef spp::sparse_hash_set<value_type> store_type;

protected:

    std::unique_ptr<store_type> _store;

public:
    
    SparseppSetStorage()
    {
        _store = std::make_unique<store_type>();
    }

    static std::shared_ptr<SparseppSetStorage> build();

    static std::shared_ptr<SparseppSetStorage> build(const typename StorageTraits<SparseppSetStorage>::params_type&);

    static std::shared_ptr<SparseppSetStorage> deserialize(std::ifstream& in);

    void serialize(std::ofstream& out);

    std::shared_ptr<SparseppSetStorage> clone() const;

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

    void save(std::string, uint16_t );

    void load(std::string, uint16_t &);

    const bool insert(value_type h);

    const count_t insert_and_query(value_type h);

    const count_t query(value_type h) const;


    byte_t ** get_raw_tables() {
        return nullptr;
    }

};


template<>
struct is_probabilistic<SparseppSetStorage> { 
    static const bool value = false;
};

template<>
struct is_counting<SparseppSetStorage> {
    static const bool value = false;
};

}
#endif
