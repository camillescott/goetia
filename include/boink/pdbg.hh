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
#include "boink/hashing/hashing_types.hh"
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
class PdBG : public kmers::KmerClient,
             public std::enable_shared_from_this<PdBG<BaseStorageType>> {
protected:

    std::unique_ptr<PartitionedStorage<BaseStorageType>>       S;
    std::shared_ptr<hashing::UKHS>                             ukhs;
    hashing::DefaultUKHSShifter                                partitioner;

public:

    typedef hashing::DefaultUKHSShifter                        shifter_type;
	typedef AssemblerMixin<PdBG<BaseStorageType>>              assembler_type;
    typedef hashing::KmerIterator<hashing::DefaultUKHSShifter> kmer_iter_type;
    typedef BaseStorageType                                    base_storage_type;

    const uint16_t partition_K;
 
    template <typename... Args>
    explicit PdBG(uint16_t  K,
                  uint16_t  partition_K,
                  std::shared_ptr<hashing::UKHS> ukhs,
                  Args&&... args)
        : KmerClient  (K),
          ukhs        (ukhs),
          partitioner (K, partition_K, ukhs),
          partition_K (partition_K)
    {
        S = std::make_unique<PartitionedStorage<BaseStorageType>>(partitioner.n_ukhs_hashes(),
                                                                  std::forward<Args>(args)...);
    }

    explicit PdBG(uint16_t K,
                  uint16_t partition_K,
                  std::shared_ptr<hashing::UKHS> ukhs,
                  std::unique_ptr<storage::PartitionedStorage<BaseStorageType>>& S)
        : KmerClient  (K),
          ukhs        (ukhs),
          partitioner (K, partition_K, ukhs),
          partition_K (partition_K),
          S(std::move(S->clone()))
    {

    }

    hashing::hash_t hash(const std::string& kmer) const {
        return partitioner.hash(kmer);
    }

    hashing::hash_t hash(const char * kmer) const {
        return partitioner.hash(kmer);
    }

    std::vector<hashing::hash_t> get_hashes(const std::string& sequence) {

        hashing::KmerIterator<hashing::DefaultUKHSShifter> iter(sequence, &partitioner);
        std::vector<hashing::hash_t> kmer_hashes;

        while(!iter.done()) {
            hashing::PartitionedHash h = iter.next();
            kmer_hashes.push_back(h.first);
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

    inline const bool insert(const std::string& kmer) {
        partitioner.set_cursor(kmer);
        partitioner.reset_unikmers();

        return S->insert(partitioner.get(), partitioner.get_front_partition());
    }

    inline const bool insert(const hashing::PartitionedHash ph) {
        return S->insert(ph.first, ph.second);
    }

    inline const storage::count_t insert_and_query(const std::string& kmer) {
        partitioner.set_cursor(kmer);
        partitioner.reset_unikmers();

        return S->insert_and_query(partitioner.get(), partitioner.get_front_partition());
    }

    inline const storage::count_t insert_and_query(const hashing::PartitionedHash ph) {
        return S->insert_and_query(ph.first, ph.second);
    }

    inline const storage::count_t query(const std::string& kmer) {
        partitioner.set_cursor(kmer);
        partitioner.reset_unikmers();

        return S->query(partitioner.get(), partitioner.get_front_partition());
    }

    inline const storage::count_t query(const hashing::PartitionedHash ph) {
        return S->query(ph.first, ph.second);
    }

    uint64_t insert_sequence(const std::string&          sequence,
                             std::vector<hashing::hash_t>&  kmer_hashes,
                             std::vector<storage::count_t>& counts) {

        hashing::KmerIterator<hashing::DefaultUKHSShifter> iter(sequence, &partitioner);

        uint64_t         n_consumed = 0;
        size_t           pos = 0;
        storage::count_t count;
        while(!iter.done()) {
            hashing::PartitionedHash h = iter.next();
            count = insert_and_query(h);

            kmer_hashes.push_back(h.first);
            counts.push_back(count);

            n_consumed += (count == 1);
            ++pos;
        }

        return n_consumed;
    }

    uint64_t insert_sequence(const std::string&      sequence,
                          std::set<hashing::hash_t>& new_kmers) {

        hashing::KmerIterator<hashing::DefaultUKHSShifter> iter(sequence, &partitioner);

        uint64_t n_consumed = 0;
        size_t pos = 0;
        bool is_new;
        while(!iter.done()) {
            hashing::PartitionedHash h = iter.next();
            if(insert(h)) {
                new_kmers.insert(h.first);
                n_consumed += is_new;
            }
            ++pos;
        }

        return n_consumed;
    }

    uint64_t insert_sequence(const std::string& sequence) {
        hashing::KmerIterator<hashing::DefaultUKHSShifter> iter(sequence, &partitioner);

        uint64_t n_consumed = 0;
        while(!iter.done()) {
            hashing::PartitionedHash h = iter.next();
            n_consumed += insert(h);
        }

        return n_consumed;
    }

    uint64_t insert_sequence_rolling(const std::string& sequence) {
        hashing::KmerIterator<hashing::DefaultUKHSShifter> iter(sequence, &partitioner);

        uint64_t n_consumed = 0;
        hashing::PartitionedHash h = iter.next();
        uint64_t          cur_pid       = h.second;
        BaseStorageType * cur_partition = S->query_partition(cur_pid);
        cur_partition->insert(h.first);

        while(!iter.done()) {
            h = iter.next();
            if (h.second != cur_pid) {
                cur_pid = h.second;
                cur_partition = S->query_partition(cur_pid);
            }
            n_consumed += cur_partition->insert(h.first);
        }

        return n_consumed;
    }

    std::vector<storage::count_t> insert_and_query_sequence(const std::string& sequence) {

        hashing::KmerIterator<hashing::DefaultUKHSShifter> iter(sequence, &partitioner);
        std::vector<storage::count_t> counts(sequence.length() - _K + 1);

        size_t pos = 0;
        while(!iter.done()) {
            hashing::PartitionedHash h = iter.next();
            counts[pos] = S->insert_and_query(h);
            ++pos;
        }

        return counts;
    }

    std::vector<storage::count_t> query_sequence(const std::string& sequence) {

        hashing::KmerIterator<hashing::DefaultUKHSShifter> iter(sequence, &partitioner);
        std::vector<storage::count_t> counts(sequence.length() - _K + 1);

        size_t pos = 0;
        while(!iter.done()) {
            hashing::PartitionedHash h = iter.next();
            counts[pos] = query(h);
            ++pos;
        }

        return counts;
    }

    std::vector<storage::count_t> query_sequence_rolling(const std::string& sequence) {

        hashing::KmerIterator<hashing::DefaultUKHSShifter> iter(sequence, &partitioner);
        std::vector<storage::count_t> counts(sequence.length() - _K + 1);
        
        hashing::PartitionedHash h      = iter.next();
        uint64_t          cur_pid       = h.second;
        BaseStorageType * cur_partition = S->query_partition(cur_pid);
        counts[0]                       = cur_partition->query(h.first);

        size_t pos = 1;
        while(!iter.done()) {
            h = iter.next();
            if (h.second != cur_pid) {
                cur_pid = h.second;
                cur_partition = S->query_partition(cur_pid);
            }
            counts[pos] = cur_partition->query(h.first);
            ++pos;
        }

        return counts;
    }

    void query_sequence(const std::string&             sequence,
                        std::vector<storage::count_t>& counts,
                        std::vector<hashing::hash_t>&  hashes) {

        hashing::KmerIterator<hashing::DefaultUKHSShifter> iter(sequence, &partitioner);

        while(!iter.done()) {
            hashing::PartitionedHash h = iter.next();
            storage::count_t result = query(h);
            counts.push_back(result);
            hashes.push_back(h.first);
        }
    }

    void query_sequence(const std::string& sequence,
                        std::vector<storage::count_t>& counts,
                        std::vector<hashing::hash_t>& hashes,
                        std::set<hashing::hash_t>& new_hashes) {

        hashing::KmerIterator<hashing::DefaultUKHSShifter> iter(sequence, &partitioner);

        while(!iter.done()) {
            hashing::PartitionedHash h = iter.next();
            auto result = query(h);
            if (result == 0) {
                new_hashes.insert(h.first);
            }
            counts.push_back(result);
            hashes.push_back(h.first);
        }
    }

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

    std::vector<hashing::kmer_t> build_left_kmers(const std::vector<hashing::shift_t>& nodes,
                                                  const std::string& root) {
        std::vector<hashing::kmer_t> kmers;
        auto _prefix = prefix(root);
        for (auto neighbor : nodes) {
            kmers.push_back(hashing::kmer_t(neighbor.hash,
                                   neighbor.symbol
                                   + _prefix));
        }
        return kmers;
    }

    std::vector<hashing::kmer_t> build_right_kmers(const std::vector<hashing::shift_t>& nodes,
                                                   const std::string& root) {
        std::vector<hashing::kmer_t> kmers;
        auto _suffix = suffix(root);
        for (auto neighbor : nodes) {
            kmers.push_back(hashing::kmer_t(neighbor.hash,
                                   _suffix
                                   + neighbor.symbol));
        }

        return kmers;
    }

    /*
    std::vector<hashing::shift_t> left_neighbors(const std::string& root) {
        auto assembler = this->get_assembler();
        assembler->set_cursor(root);
        return assembler->filter_nodes(assembler->gather_left());
    }

    std::vector<hashing::kmer_t> left_neighbor_kmers(const std::string& root) {
        auto filtered = left_neighbors(root);
        return build_left_kmers(filtered, root);
    }

    std::vector<hashing::shift_t> right_neighbors(const std::string& root) {
        auto assembler = this->get_assembler();
        assembler->set_cursor(root);
        return assembler->filter_nodes(assembler->gather_right());
    }

    std::vector<hashing::kmer_t> right_neighbor_kmers(const std::string& root) {
        auto filtered = right_neighbors(root);
        return build_right_kmers(filtered, root);
    }

    std::pair<std::vector<hashing::shift_t>,
              std::vector<hashing::shift_t>> neighbors(const std::string& root) {
        auto assembler = this->get_assembler();
        assembler->set_cursor(root);
        auto lfiltered = assembler->filter_nodes(assembler->gather_left());
        auto rfiltered = assembler->filter_nodes(assembler->gather_right());
        return std::make_pair(lfiltered, rfiltered);
    }

    std::pair<std::vector<hashing::kmer_t>,
              std::vector<hashing::kmer_t>> neighbor_kmers(const std::string& root) {
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

    std::shared_ptr<hashing::KmerIterator<hashing::DefaultUKHSShifter>> 
    get_hash_iter(const std::string& sequence) {
        return std::make_shared<hashing::KmerIterator<hashing::DefaultUKHSShifter>>(sequence, &partitioner);
    }

    /*
    std::shared_ptr<assembler_type> get_assembler() {
        auto ptr = this->shared_from_this();
        return std::make_shared<assembler_type>(ptr, partitioner);
    }
    */

    std::vector<size_t> get_partition_counts() {
        return S->get_partition_counts();
    }
};

typedef PdBG<storage::SparseppSetStorage> DefaultPdBG;

}
#endif
