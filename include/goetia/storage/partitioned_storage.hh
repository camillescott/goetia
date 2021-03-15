/**
 * (c) Camille Scott, 2019
 * File   : partitioned_storage.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */
/* partitioned_storage.hh -- storage classes for the goetia dbg
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef GOETIA_PARTITIONEDSTORAGE_HH
#define GOETIA_PARTITIONEDSTORAGE_HH

#include "goetia/goetia.hh"
#include "goetia/storage/storage.hh"
#include "goetia/storage/storage_types.hh"
#include "sparsepp/spp.h"

#include <vector>

namespace goetia {
namespace storage {


template <class BaseStorageType>
class PartitionedStorage : public Storage<uint64_t> {

protected:
    std::vector<std::shared_ptr<BaseStorageType>> partitions;
    const uint64_t                                n_partitions;

public:

    typedef uint64_t        value_type;
    typedef BaseStorageType base_storage_type;

    PartitionedStorage (const uint64_t n_partitions)
        : PartitionedStorage(n_partitions, StorageTraits<BaseStorageType>::default_params)
    {
    }
    
    PartitionedStorage (const uint64_t n_partitions,
                        const typename StorageTraits<BaseStorageType>::params_type& params)
        : n_partitions(n_partitions)
    {
        for (size_t i = 0; i < n_partitions; ++i) {
            partitions.push_back(
                base_storage_type::build(params)
            );
        }
    }
    
    template <typename... Args>
    PartitionedStorage(const uint64_t  n_partitions,
                       Args&&...       args)
        : n_partitions(n_partitions)
    {
        for (size_t i = 0; i < n_partitions; ++i) {
            partitions.push_back(
                std::make_shared<base_storage_type>(std::forward<Args>(args)...)
            );
        }
    }

    PartitionedStorage(const uint64_t n_partitions,
                       base_storage_type* S)
        : n_partitions(n_partitions)
    {
        for (size_t i = 0; i < n_partitions; ++i) {
            partitions.push_back(std::move(S->clone()));
        }
    }

    std::shared_ptr<PartitionedStorage<BaseStorageType>> clone() const {
        return std::make_shared<PartitionedStorage<BaseStorageType>>(n_partitions,
                                                                     partitions.front().get());

    }

    void reset() {
        for (size_t i = 0; i < n_partitions; ++i) {
            partitions[i]->reset();
        }
    }

    //std::vector<uint64_t> get_tablesizes() const {
    //    return partitions.front()->get_tablesizes();
    //}

    const uint64_t n_unique_kmers() const {
        uint64_t sum = 0;
        for (auto& partition : partitions) {
            sum += partition->n_unique_kmers();
        }
        return sum;
    }

    //const uint64_t n_tables() const {
    //    return partitions.front()->n_tables();
   // }

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

    inline const bool insert(value_type h, uint64_t partition) {
        return query_partition(partition)->insert(h);
    }

    inline const count_t insert_and_query(value_type h, uint64_t partition) {
        auto partition_store = query_partition(partition);
        return partition_store->insert_and_query(h);
    }

    inline const count_t query(value_type h, uint64_t partition) {
        return query_partition(partition)->query(h);
    }


    BaseStorageType * query_partition(uint64_t partition) {
        if (partition < n_partitions) {
            return partitions[partition].get();
        } else {
            throw GoetiaException("Invalid storage partition: " + std::to_string(partition));
        }
    }

    byte_t ** get_raw_tables() {
        return nullptr;
    }


    const bool insert(value_type khash ) {
        throw GoetiaException("Method not available!");
    }

    const storage::count_t insert_and_query(value_type khash) {
        throw GoetiaException("Method not available!");
    }

    const storage::count_t query(value_type khash) const {
        throw GoetiaException("Method not available!");
    }

    std::vector<size_t> get_partition_counts() {
        std::vector<size_t> counts;
        for (auto& partition : partitions) {
            counts.push_back(partition->n_unique_kmers());
        }
        return counts;
    }

    void * get_partition_counts_as_buffer() {
        size_t * counts = (size_t *) malloc(sizeof(size_t) * n_partitions);
        for (size_t pidx = 0; pidx < n_partitions; ++pidx) {
            counts[pidx] = partitions[pidx]->n_unique_kmers();
        }
        return counts;
    }
};


template<class StorageType>
struct is_probabilistic<PartitionedStorage<StorageType>> { 
      static const bool value = is_probabilistic<StorageType>::value;
};


extern template class goetia::storage::PartitionedStorage<goetia::storage::BitStorage>;
extern template class goetia::storage::PartitionedStorage<goetia::storage::ByteStorage>;
extern template class goetia::storage::PartitionedStorage<goetia::storage::NibbleStorage>;
extern template class goetia::storage::PartitionedStorage<goetia::storage::QFStorage>;
extern template class goetia::storage::PartitionedStorage<goetia::storage::SparseppSetStorage>;
extern template class goetia::storage::PartitionedStorage<goetia::storage::PHMapStorage>;
extern template class goetia::storage::PartitionedStorage<goetia::storage::BTreeStorage>;

}
}
#endif
