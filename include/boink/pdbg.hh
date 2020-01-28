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
#include "boink/meta.hh"
#include "boink/traversal.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/hashextender.hh"
#include "boink/hashing/canonical.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/hashing/shifter_types.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/storage_types.hh"
#include "boink/storage/partitioned_storage.hh"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace boink {

using storage::PartitionedStorage;

template <class BaseStorageType,
          class ShifterType = hashing::FwdUnikmerShifter>
class PdBG {

public:
    
    typedef PdBG<BaseStorageType, ShifterType>                 graph_type;
    typedef ShifterType                                        shifter_type;
    typedef dBGWalker<PdBG<BaseStorageType, ShifterType>>      walker_type;
    typedef typename shifter_type::ukhs_type                   ukhs_type;
    _boink_model_typedefs_from_shiftertype(shifter_type)
    _boink_walker_typedefs_from_graphtype(walker_type)

    typedef BaseStorageType                                    base_storage_type;

protected:

    std::shared_ptr<PartitionedStorage<BaseStorageType>> S;
    std::shared_ptr<ukhs_type>                           ukhs;
    extender_type                                        partitioner;

public:

    const uint16_t K;
    const uint16_t partition_K;
 
    template <typename... Args>
    explicit PdBG(uint16_t  K,
                  uint16_t  partition_K,
                  std::shared_ptr<ukhs_type>& ukhs,
                  Args&&... args)
        : K           (K),
          ukhs        (ukhs),
          partitioner (K, partition_K, ukhs),
          partition_K (partition_K)
    {
        S = std::make_shared<PartitionedStorage<BaseStorageType>>(ukhs->n_hashes(),
                                                                  std::forward<Args>(args)...);
    }

    explicit PdBG(uint16_t K,
                  uint16_t partition_K,
                  std::shared_ptr<ukhs_type>& ukhs,
                  std::shared_ptr<storage::PartitionedStorage<BaseStorageType>> S)
        : K           (K),
          ukhs        (ukhs),
          partitioner (K, partition_K, ukhs),
          partition_K (partition_K),
          S(S->clone())
    {

    }

    inline hash_type hash(const std::string& kmer) const {
        return shifter_type::hash(kmer, K, partition_K, ukhs);
    }

    inline hash_type hash(const char * kmer) const {
        return shifter_type::hash(kmer, K, partition_K, ukhs);
    }

    inline std::vector<hash_type> get_hashes(const std::string& sequence) {

        hashing::KmerIterator<extender_type> iter(sequence, &partitioner);
        std::vector<hash_type> kmer_hashes;

        while(!iter.done()) {
            hash_type h = iter.next();
            kmer_hashes.push_back(h);
        }

        return kmer_hashes;
    }

    inline std::shared_ptr<graph_type> clone() {
        return std::make_shared<graph_type>(K,
                                            partition_K,
                                            ukhs,
                                            S);
    }

    inline const bool insert(const std::string& kmer) {

        partitioner.set_cursor(kmer);
        auto bh = partitioner.get();

        return S->insert(bh.value(), bh.minimizer.partition);
    }

    inline const bool insert(const hash_type& h) {
        return S->insert(h.value(), h.minimizer.partition);
    }

    inline const storage::count_t insert_and_query(const std::string& kmer) {

        partitioner.set_cursor(kmer);
        auto h = partitioner.get();

        return S->insert_and_query(h.value(), h.minimizer.partition);
    }

    inline const storage::count_t insert_and_query(const hash_type& h) {
        return S->insert_and_query(h.value(), h.minimizer.partition);
    }

    inline const storage::count_t query(const std::string& kmer) {
        auto h = hash(kmer);
        return S->query(h.value(), h.minimizer.partition);
    }

    inline const storage::count_t query(const hash_type& h) {
        return S->query(h.value(), h.minimizer.partition);
    }

    inline const uint64_t insert_sequence(const std::string&             sequence,
                                   std::vector<hash_type>&        kmer_hashes,
                                   std::vector<storage::count_t>& counts) {

        hashing::KmerIterator<extender_type> iter(sequence, &partitioner);

        uint64_t         n_consumed = 0;
        size_t           pos = 0;
        storage::count_t count;
        while(!iter.done()) {
            auto h = iter.next();
            count = insert_and_query(h);

            kmer_hashes.push_back(h);
            counts.push_back(count);

            n_consumed += (count == 1);
            ++pos;
        }

        return n_consumed;
    }

    inline const uint64_t insert_sequence(const std::string&   sequence,
                                   std::set<hash_type>& new_kmers) {

        hashing::KmerIterator<extender_type> iter(sequence, &partitioner);

        uint64_t n_consumed = 0;
        size_t pos = 0;
        bool is_new;
        while(!iter.done()) {
            auto h = iter.next();
            if(insert(h)) {
                new_kmers.insert(h);
                n_consumed += is_new;
            }
            ++pos;
        }

        return n_consumed;
    }

    inline const uint64_t insert_sequence(const std::string& sequence) {
        hashing::KmerIterator<extender_type> iter(sequence, &partitioner);

        uint64_t n_consumed = 0;
        while(!iter.done()) {
            hash_type h = iter.next();
            n_consumed += insert(h);
        }

        return n_consumed;
    }

    inline const uint64_t insert_sequence_rolling(const std::string& sequence) {
        hashing::KmerIterator<extender_type> iter(sequence, &partitioner);

        uint64_t          n_consumed    = 0;
        auto              h             = iter.next();
        uint64_t          cur_pid       = h.minimizer.partition;
        BaseStorageType * cur_partition = S->query_partition(cur_pid);
        cur_partition->insert(h.value());

        while(!iter.done()) {
            h = iter.next();
            if (h.minimizer.partition != cur_pid) {
                cur_pid = h.minimizer.partition;
                cur_partition = S->query_partition(cur_pid);
            }
            n_consumed += cur_partition->insert(h.value());
        }

        return n_consumed;
    }


    inline std::vector<storage::count_t> insert_and_query_sequence(const std::string& sequence) {

        hashing::KmerIterator<extender_type> iter(sequence, &partitioner);
        std::vector<storage::count_t> counts(sequence.length() - K + 1);

        size_t pos = 0;
        while(!iter.done()) {
            auto h = iter.next();
            counts[pos] = insert_and_query(h);
            ++pos;
        }

        return counts;
    }


    inline std::vector<storage::count_t> query_sequence(const std::string& sequence) {

        hashing::KmerIterator<extender_type> iter(sequence, &partitioner);
        std::vector<storage::count_t> counts(sequence.length() - K + 1);

        size_t pos = 0;
        while(!iter.done()) {
            auto h = iter.next();
            counts[pos] = query(h);
            ++pos;
        }

        return counts;
    }


    inline std::vector<storage::count_t> query_sequence_rolling(const std::string& sequence) {

        hashing::KmerIterator<extender_type> iter(sequence, &partitioner);
        std::vector<storage::count_t> counts(sequence.length() - K + 1);
        
        hash_type         h             = iter.next();
        uint64_t          cur_pid       = h.minimizer.partition;
        BaseStorageType * cur_partition = S->query_partition(cur_pid);
        counts[0]                       = cur_partition->query(h.value());

        size_t pos = 1;
        while(!iter.done()) {
            h = iter.next();
            if (h.minimizer.partition != cur_pid) {
                cur_pid = h.minimizer.partition;
                cur_partition = S->query_partition(cur_pid);
            }
            counts[pos] = cur_partition->query(h.value());
            ++pos;
        }

        return counts;
    }


    inline void query_sequence(const std::string&             sequence,
                        std::vector<storage::count_t>& counts,
                        std::vector<hash_type>&        hashes) {

        hashing::KmerIterator<extender_type> iter(sequence, &partitioner);

        while(!iter.done()) {
            auto h = iter.next();
            storage::count_t result = query(h);
            counts.push_back(result);
            hashes.push_back(h);
        }
    }


    inline void query_sequence(const std::string&             sequence,
                        std::vector<storage::count_t>& counts,
                        std::vector<hash_type>&        hashes,
                        std::set<hash_type>&           new_hashes) {

        hashing::KmerIterator<extender_type> iter(sequence, &partitioner);

        while(!iter.done()) {
            auto h = iter.next();
            auto result = query(h);
            if (result == 0) {
                new_hashes.insert(h);
            }
            counts.push_back(result);
            hashes.push_back(h);
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
        return kmer.substr(kmer.length() - this->K + 1);
    }

    const std::string prefix(const std::string& kmer) {
        return kmer.substr(0, this->K - 1);
    }

    std::vector<kmer_type> build_left_kmers(const std::vector<shift_type<hashing::DIR_LEFT>>& nodes,
                                            const std::string& root) {
        std::vector<kmer_type> kmers;
        auto _prefix = prefix(root);
        for (auto neighbor : nodes) {
            kmers.push_back(kmer_type(neighbor,
                                      neighbor.symbol + _prefix));
        }
        return kmers;
    }

    std::vector<kmer_type> build_right_kmers(const std::vector<shift_type<hashing::DIR_RIGHT>>& nodes,
                                             const std::string& root) {
        std::vector<kmer_type> kmers;
        auto _suffix = suffix(root);
        for (auto neighbor : nodes) {
            kmers.push_back(kmer_type(neighbor,
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
        S->save(filename, K);
    }

    void load(std::string filename) {
        uint16_t ksize = K;
        S->load(filename, ksize);
    }

    void reset() {
        S->reset();
    }

    std::vector<size_t> get_partition_counts() {
        return S->get_partition_counts();
    }
};


extern template class PdBG<storage::BitStorage, hashing::FwdUnikmerShifter>;
extern template class PdBG<storage::BitStorage, hashing::CanUnikmerShifter>;

extern template class PdBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>;
extern template class PdBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>;

extern template class PdBG<storage::ByteStorage, hashing::FwdUnikmerShifter>;
extern template class PdBG<storage::ByteStorage, hashing::CanUnikmerShifter>;

extern template class PdBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>;
extern template class PdBG<storage::NibbleStorage, hashing::CanUnikmerShifter>;

extern template class PdBG<storage::QFStorage, hashing::FwdUnikmerShifter>;
extern template class PdBG<storage::QFStorage, hashing::CanUnikmerShifter>;

}
#endif
