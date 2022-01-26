/**
 * (c) Camille Scott, 2021
 * File   : phmapstorage.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 12.03.2021
 */

#ifndef GOETIA_PHMAPSTORAGE_HH
#define GOETIA_PHMAPSTORAGE_HH

#include "goetia/meta.hh"
#include "goetia/storage/storage.hh"
#include "goetia/storage/phmap/phmap.h"

#include <cstdint>
#include <fstream>
#include <memory>
#include <tuple>
#include <vector>

namespace goetia {

class PHMapStorage;

template<>
struct StorageTraits<PHMapStorage> {
    static constexpr bool is_probabilistic = false;
    static constexpr bool is_counting      = false;

    typedef std::tuple<bool> params_type;
    static constexpr params_type default_params = std::make_tuple(0);
};


class PHMapStorage : public Storage<uint64_t>,
                           public Tagged<PHMapStorage> {

public:
    
    using Storage<uint64_t>::value_type;
    using Traits = StorageTraits<PHMapStorage>;
    typedef phmap::parallel_flat_hash_set<value_type> store_type;

protected:

    std::unique_ptr<store_type> _store;

public:
    
    PHMapStorage()
    {
        _store = std::make_unique<store_type>();
    }

    static std::shared_ptr<PHMapStorage> build();

    static std::shared_ptr<PHMapStorage> build(const typename StorageTraits<PHMapStorage>::params_type&);

    static std::shared_ptr<PHMapStorage> deserialize(std::ifstream& in);

    void serialize(std::ofstream& out);

    std::shared_ptr<PHMapStorage> clone() const;

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

    const inline bool insert(value_type h) {
        auto result = _store->insert(h);
        // the second in the returned pair reports that the insert
        // took place ie the hash was new
        return result.second;
    }

    const count_t insert_and_query(value_type h);

    const count_t query(value_type h) const;


    byte_t ** get_raw_tables() {
        return nullptr;
    }

};


}
#endif
