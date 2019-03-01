/* partitioned_storage.hh -- storage classes for the boink dbg
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_PARTITIONEDSTORAGE_HH
#define BOINK_PARTITIONEDSTORAGE_HH

#include "boink/boink.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/storage/storage.hh"
#include "sparsepp/spp.h"

#include <vector>

namespace boink {
namespace storage {


template <class BaseStorageType>
class PartitionedStorage : public Storage {

protected:
    std::vector<std::unique_ptr<BaseStorageType>> partitions;
    const uint64_t                                n_partitions;

public:

    typedef BaseStorageType base_storage_type;
    
    template <typename... Args>
    PartitionedStorage(const uint64_t  n_partitions,
                       Args&&...       args)
        : n_partitions(n_partitions)
    {
        for (size_t i = 0; i < n_partitions; ++i) {
            partitions.push_back(
                std::make_unique<BaseStorageType>(std::forward<Args>(args)...)
            );
        }
    }

    PartitionedStorage(const uint64_t n_partitions,
                       BaseStorageType* S)
        : n_partitions(n_partitions)
    {
        for (size_t i = 0; i < n_partitions; ++i) {
            partitions.push_back(std::move(S->clone()));
        }
    }

    std::unique_ptr<PartitionedStorage<BaseStorageType>> clone() const {
        return std::make_unique<PartitionedStorage<BaseStorageType>>(n_partitions,
                                                                     partitions.front().get());

    }

    void reset() {
        for (size_t i = 0; i < n_partitions; ++i) {
            partitions[i]->reset();
        }
    }

    std::vector<uint64_t> get_tablesizes() const {
        return partitions.front()->get_tablesizes();
    }

    const uint64_t n_unique_kmers() const {
        uint64_t sum = 0;
        for (auto& partition : partitions) {
            sum += partition->n_unique_kmers();
        }
        return sum;
    }

    const uint64_t n_tables() const {
        return partitions.front()->n_tables();
    }

    const uint64_t n_occupied() const {
        return partitions.front()->n_occupied();
    }

    const uint64_t n_partition_stores() const {
        return n_partitions;
    }

    template<typename Dummy = double>
    auto estimated_fp() 
    -> std::enable_if_t<storage::is_probabilistic<BaseStorageType>::value, Dummy>
    {
        double sum = 0;
        for (auto& partition : partitions) {
            sum += partition->estimated_fp();
        }
        return sum / (double)n_partition_stores();
    }

    void save(std::string, uint16_t ) {
     
    }

    void load(std::string, uint16_t &) {

    }

    inline const bool insert(hashing::hash_t h, uint64_t partition) {
        return query_partition(partition)->insert(h);
    }

    inline const count_t insert_and_query(hashing::hash_t h, uint64_t partition) {
        auto partition_store = query_partition(partition);
        return partition_store->insert_and_query(h);
    }

    inline const count_t query(hashing::hash_t h, uint64_t partition) {
        return query_partition(partition)->query(h);
    }


    BaseStorageType * query_partition(uint64_t partition) {
        if (partition < n_partitions) {
            return partitions[partition].get();
        } else {
            throw BoinkException("Invalid storage partition: " + std::to_string(partition));
        }
    }

    byte_t ** get_raw_tables() {
        return nullptr;
    }


    const bool insert(hashing::hash_t khash ) {
        throw BoinkException("Method not available!");
    }

    const storage::count_t insert_and_query(hashing::hash_t khash) {
        throw BoinkException("Method not available!");
    }

    const storage::count_t query(hashing::hash_t khash) const {
        throw BoinkException("Method not available!");
    }

    std::vector<size_t> get_partition_counts() {
        std::vector<size_t> counts;
        for (auto& partition : partitions) {
            counts.push_back(partition->n_unique_kmers());
        }
        return counts;
    }
};


template<>
template<class StorageType>
struct is_probabilistic<PartitionedStorage<StorageType>> { 
      static const bool value = is_probabilistic<StorageType>::value;
};

}
}
#endif
