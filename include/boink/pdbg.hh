/**
 * (c) Camille Scott, 2019
 * File   : pdbg.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */
/* dbg.hh -- Generic de Bruijn graph with partitioned storagee
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_PDBG_HH
#define BOINK_PDBG_HH

#include "boink/boink.hh"
#include "boink/assembly.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/hashshifter.hh"
#include "boink/kmers/kmerclient.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/sparseppstorage.hh"
#include "boink/storage/partitioned_storage.hh"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace boink {

using storage::PartitionedStorage;

template <class BaseStorageType>
class PdBG : public kmers::KmerClient {
protected:

    std::shared_ptr<PartitionedStorage<BaseStorageType>>       S;
    std::shared_ptr<hashing::UKHS::Map>                        ukhs;
    hashing::UKHS::LazyShifter                                 partitioner;

public:

    typedef typename hashing::UKHS::LazyShifter                          shifter_type;
    typedef typename shifter_type::hash_type                             hash_type;
    typedef typename shifter_type::shift_type                            shift_type;
    typedef typename shifter_type::kmer_type                             kmer_type;

	typedef Traverse<PdBG<BaseStorageType>>                              traversal_type;
    typedef hashing::KmerIterator<typename hashing::UKHS::LazyShifter>   kmer_iter_type;
    typedef BaseStorageType                                              base_storage_type;

    const uint16_t partition_K;
 
    template <typename... Args>
    explicit PdBG(uint16_t  K,
                  uint16_t  partition_K,
                  std::shared_ptr<hashing::UKHS::Map> ukhs,
                  Args&&... args)
        : KmerClient  (K),
          ukhs        (ukhs),
          partitioner (K, partition_K, ukhs),
          partition_K (partition_K)
    {
        S = std::make_shared<PartitionedStorage<BaseStorageType>>(partitioner.n_ukhs_hashes(),
                                                                  std::forward<Args>(args)...);
    }

    explicit PdBG(uint16_t K,
                  uint16_t partition_K,
                  std::shared_ptr<hashing::UKHS::Map> ukhs,
                  std::shared_ptr<storage::PartitionedStorage<BaseStorageType>> S)
        : KmerClient  (K),
          ukhs        (ukhs),
          partitioner (K, partition_K, ukhs),
          partition_K (partition_K),
          S(S->clone())
    {

    }

    hash_type hash(const std::string& kmer) const {
        return partitioner.hash(kmer);
    }

    hash_type hash(const char * kmer) const {
        return partitioner.hash(kmer);
    }

    std::vector<hash_type> get_hashes(const std::string& sequence) {

        kmer_iter_type iter(sequence, partitioner);
        std::vector<hash_type> kmer_hashes;

        while(!iter.done()) {
            hash_type h = iter.next();
            kmer_hashes.push_back(h);
        }

        return kmer_hashes;
    }

    std::shared_ptr<PdBG<BaseStorageType>> clone() {
        return std::make_shared<PdBG<BaseStorageType>>(
                   this->_K,
                   partition_K,
                   ukhs,
                   S
               );
    }

    const bool insert(const std::string& kmer);

    const bool insert(const hash_type& h);

    const storage::count_t insert_and_query(const std::string& kmer);

    const storage::count_t insert_and_query(const hash_type& h);

    const storage::count_t query(const std::string& kmer);

    const storage::count_t query(const hash_type& h);

    uint64_t insert_sequence(const std::string&             sequence,
                             std::vector<hash_type>&        kmer_hashes,
                             std::vector<storage::count_t>& counts);

    uint64_t insert_sequence(const std::string&   sequence,
                             std::set<hash_type>& new_kmers);

    uint64_t insert_sequence(const std::string& sequence);

    uint64_t insert_sequence_rolling(const std::string& sequence);

    std::vector<storage::count_t> insert_and_query_sequence(const std::string& sequence);

    std::vector<storage::count_t> query_sequence(const std::string& sequence);

    std::vector<storage::count_t> query_sequence_rolling(const std::string& sequence);

    void query_sequence(const std::string&             sequence,
                        std::vector<storage::count_t>& counts,
                        std::vector<hash_type>&         hashes);

    void query_sequence(const std::string&             sequence,
                        std::vector<storage::count_t>& counts,
                        std::vector<hash_type>&        hashes,
                        std::set<hash_type>&           new_hashes);

    uint64_t n_unique() const {
        return S->n_unique_kmers();
    }

    uint64_t n_occupied() const {
        return S->n_occupied();
    }

    uint64_t n_partitions() const {
        return S->n_partition_stores();
    }

    const std::string suffix(const std::string& kmer) {
        return kmer.substr(kmer.length() - this->_K + 1);
    }

    const std::string prefix(const std::string& kmer) {
        return kmer.substr(0, this->_K - 1);
    }

    std::vector<kmer_type> build_left_kmers(const std::vector<shift_type>& nodes,
                                            const std::string& root) {
        std::vector<kmer_type> kmers;
        auto _prefix = prefix(root);
        for (auto neighbor : nodes) {
            kmers.push_back(kmer_type(neighbor.hash,
                                      neighbor.symbol + _prefix));
        }
        return kmers;
    }

    std::vector<kmer_type> build_right_kmers(const std::vector<shift_type>& nodes,
                                             const std::string& root) {
        std::vector<kmer_type> kmers;
        auto _suffix = suffix(root);
        for (auto neighbor : nodes) {
            kmers.push_back(kmer_type(neighbor.hash,
                                      _suffix + neighbor.symbol));
        }

        return kmers;
    }

    /*
    std::vector<shift_type> left_neighbors(const std::string& root) {
        auto assembler = this->get_assembler();
        assembler->set_cursor(root);
        return assembler->filter_nodes(assembler->gather_left());
    }

    std::vector<kmer_type> left_neighbor_kmers(const std::string& root) {
        auto filtered = left_neighbors(root);
        return build_left_kmers(filtered, root);
    }

    std::vector<shift_type> right_neighbors(const std::string& root) {
        auto assembler = this->get_assembler();
        assembler->set_cursor(root);
        return assembler->filter_nodes(assembler->gather_right());
    }

    std::vector<kmer_type> right_neighbor_kmers(const std::string& root) {
        auto filtered = right_neighbors(root);
        return build_right_kmers(filtered, root);
    }

    std::pair<std::vector<shift_type>,
              std::vector<shift_type>> neighbors(const std::string& root) {
        auto assembler = this->get_assembler();
        assembler->set_cursor(root);
        auto lfiltered = assembler->filter_nodes(assembler->gather_left());
        auto rfiltered = assembler->filter_nodes(assembler->gather_right());
        return std::make_pair(lfiltered, rfiltered);
    }

    std::pair<std::vector<kmer_type>,
              std::vector<kmer_type>> neighbor_kmers(const std::string& root) {
        auto filtered = neighbors(root);
        return std::make_pair(build_left_kmers(filtered.first, root),
                              build_right_kmers(filtered.second, root));
    }
    */

    void save(std::string filename) {
        S->save(filename, _K);
    }

    void load(std::string filename) {
        uint16_t ksize = _K;
        S->load(filename, ksize);
    }

    void reset() {
        S->reset();
    }

    std::shared_ptr<hashing::KmerIterator<hashing::UKHS::LazyShifter>> get_hash_iter(const std::string& sequence) {
        return std::make_shared<hashing::KmerIterator<hashing::UKHS::LazyShifter>>(sequence, &partitioner);
    }

    std::vector<size_t> get_partition_counts() {
        return S->get_partition_counts();
    }
};

typedef PdBG<storage::SparseppSetStorage> DefaultPdBG;

}
#endif
