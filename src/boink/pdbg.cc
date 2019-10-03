/**
 * (c) Camille Scott, 2019
 * File   : pdbg.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */
#include "boink/pdbg.hh"

#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/sparseppstorage.hh"

#include <memory>

namespace boink {


template <class BaseStorageType>
const bool 
PdBG<BaseStorageType>::insert(const std::string& kmer) {

    partitioner.set_cursor(kmer);
    auto bh = partitioner.get();

    return S->insert(bh.hash, bh.unikmer.partition);
}


template <class BaseStorageType>
const bool 
PdBG<BaseStorageType>::insert(const hash_type& h) {
    return S->insert(h.hash, h.unikmer.partition);
}


template <class BaseStorageType>
const storage::count_t 
PdBG<BaseStorageType>::insert_and_query(const std::string& kmer) {

    partitioner.set_cursor(kmer);
    auto h = partitioner.get();

    return S->insert_and_query(h.hash, h.unikmer.partition);
}


template <class BaseStorageType>
const storage::count_t 
PdBG<BaseStorageType>::insert_and_query(const hash_type& h) {
    return S->insert_and_query(h.hash, h.unikmer.partition);
}


template <class BaseStorageType>
const storage::count_t 
PdBG<BaseStorageType>::query(const std::string& kmer) {

    partitioner.set_cursor(kmer);
    auto h = partitioner.get();

    return S->query(h.hash, h.unikmer.partition);
}


template <class BaseStorageType>
const storage::count_t 
PdBG<BaseStorageType>::query(const hash_type& h) {
    return S->query(h.hash, h.unikmer.partition);
}


template <class BaseStorageType>
uint64_t 
PdBG<BaseStorageType>::insert_sequence(const std::string&             sequence,
                                       std::vector<hash_type>&        kmer_hashes,
                                       std::vector<storage::count_t>& counts) {

    kmer_iter_type iter(sequence, &partitioner);

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


template <class BaseStorageType>
uint64_t 
PdBG<BaseStorageType>::insert_sequence(const std::string&   sequence,
                                       std::set<hash_type>& new_kmers) {

    kmer_iter_type iter(sequence, &partitioner);

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


template <class BaseStorageType>
uint64_t 
PdBG<BaseStorageType>::insert_sequence(const std::string& sequence) {
    kmer_iter_type iter(sequence, &partitioner);

    uint64_t n_consumed = 0;
    while(!iter.done()) {
        auto h = iter.next();
        n_consumed += insert(h);
    }

    return n_consumed;
}


template <class BaseStorageType>
uint64_t 
PdBG<BaseStorageType>::insert_sequence_rolling(const std::string& sequence) {
    kmer_iter_type iter(sequence, &partitioner);

    uint64_t          n_consumed    = 0;
    auto              h             = iter.next();
    uint64_t          cur_pid       = h.unikmer.partition;
    BaseStorageType * cur_partition = S->query_partition(cur_pid);
    cur_partition->insert(h.hash);

    while(!iter.done()) {
        h = iter.next();
        if (h.unikmer.partition != cur_pid) {
            cur_pid = h.unikmer.partition;
            cur_partition = S->query_partition(cur_pid);
        }
        n_consumed += cur_partition->insert(h.hash);
    }

    return n_consumed;
}


template <class BaseStorageType>
std::vector<storage::count_t> 
PdBG<BaseStorageType>::insert_and_query_sequence(const std::string& sequence) {

    kmer_iter_type iter(sequence, &partitioner);
    std::vector<storage::count_t> counts(sequence.length() - _K + 1);

    size_t pos = 0;
    while(!iter.done()) {
        auto h = iter.next();
        counts[pos] = insert_and_query(h);
        ++pos;
    }

    return counts;
}


template <class BaseStorageType>
std::vector<storage::count_t> 
PdBG<BaseStorageType>::query_sequence(const std::string& sequence) {

    kmer_iter_type iter(sequence, &partitioner);
    std::vector<storage::count_t> counts(sequence.length() - _K + 1);

    size_t pos = 0;
    while(!iter.done()) {
        auto h = iter.next();
        counts[pos] = query(h);
        ++pos;
    }

    return counts;
}


template <class BaseStorageType>
std::vector<storage::count_t> 
PdBG<BaseStorageType>::query_sequence_rolling(const std::string& sequence) {

    kmer_iter_type iter(sequence, &partitioner);
    std::vector<storage::count_t> counts(sequence.length() - _K + 1);
    
    hash_type         h             = iter.next();
    uint64_t          cur_pid       = h.unikmer.partition;
    BaseStorageType * cur_partition = S->query_partition(cur_pid);
    counts[0]                       = cur_partition->query(h.hash);

    size_t pos = 1;
    while(!iter.done()) {
        h = iter.next();
        if (h.unikmer.partition != cur_pid) {
            cur_pid = h.unikmer.partition;
            cur_partition = S->query_partition(cur_pid);
        }
        counts[pos] = cur_partition->query(h.hash);
        ++pos;
    }

    return counts;
}


template <class BaseStorageType>
void 
PdBG<BaseStorageType>::query_sequence(const std::string&             sequence,
                                      std::vector<storage::count_t>& counts,
                                      std::vector<hash_type>&        hashes) {

    kmer_iter_type iter(sequence, &partitioner);

    while(!iter.done()) {
        auto h = iter.next();
        storage::count_t result = query(h);
        counts.push_back(result);
        hashes.push_back(h);
    }
}


template <class BaseStorageType>
void 
PdBG<BaseStorageType>::query_sequence(const std::string&             sequence,
                                      std::vector<storage::count_t>& counts,
                                      std::vector<hash_type>&        hashes,
                                      std::set<hash_type>&           new_hashes) {

    kmer_iter_type iter(sequence, &partitioner);

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

}

template class boink::PdBG<boink::storage::BitStorage<>>;
template class boink::PdBG<boink::storage::ByteStorage>;
template class boink::PdBG<boink::storage::NibbleStorage<>>;
template class boink::PdBG<boink::storage::QFStorage<>>;
template class boink::PdBG<boink::storage::SparseppSetStorage<>>;

