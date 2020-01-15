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


template <class BaseStorageType,
          class BaseShifterType>
const bool 
PdBG<BaseStorageType, BaseShifterType>
::insert(const std::string& kmer) {

    partitioner.set_cursor(kmer);
    auto bh = partitioner.get();

    return S->insert(bh.value(), bh.minimizer.partition);
}


template <class BaseStorageType,
          class BaseShifterType>
const bool 
PdBG<BaseStorageType, BaseShifterType>
::insert(const hash_type& h) {
    return S->insert(h.value(), h.minimizer.partition);
}


template <class BaseStorageType,
          class BaseShifterType>
const storage::count_t 
PdBG<BaseStorageType, BaseShifterType>
::insert_and_query(const std::string& kmer) {

    partitioner.set_cursor(kmer);
    auto h = partitioner.get();

    return S->insert_and_query(h.value(), h.minimizer.partition);
}


template <class BaseStorageType,
          class BaseShifterType>
const storage::count_t 
PdBG<BaseStorageType, BaseShifterType>
::insert_and_query(const hash_type& h) {
    return S->insert_and_query(h.value(), h.minimizer.partition);
}


template <class BaseStorageType,
          class BaseShifterType>
const storage::count_t 
PdBG<BaseStorageType, BaseShifterType>
::query(const std::string& kmer) {

    partitioner.set_cursor(kmer);
    auto h = partitioner.get();

    return S->query(h.value(), h.minimizer.partition);
}


template <class BaseStorageType,
          class BaseShifterType>
const storage::count_t 
PdBG<BaseStorageType, BaseShifterType>
::query(const hash_type& h) {
    return S->query(h.value(), h.minimizer.partition);
}


template <class BaseStorageType,
          class BaseShifterType>
uint64_t 
PdBG<BaseStorageType, BaseShifterType>
::insert_sequence(const std::string&             sequence,
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


template <class BaseStorageType,
          class BaseShifterType>
uint64_t 
PdBG<BaseStorageType, BaseShifterType>
::insert_sequence(const std::string&   sequence,
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


template <class BaseStorageType,
          class BaseShifterType>
uint64_t 
PdBG<BaseStorageType, BaseShifterType>
::insert_sequence(const std::string& sequence) {
    hashing::KmerIterator<extender_type> iter(sequence, &partitioner);

    uint64_t n_consumed = 0;
    while(!iter.done()) {
        hash_type h = iter.next();
        n_consumed += insert(h);
    }

    return n_consumed;
}


template <class BaseStorageType,
          class BaseShifterType>
uint64_t 
PdBG<BaseStorageType, BaseShifterType>
::insert_sequence_rolling(const std::string& sequence) {
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


template <class BaseStorageType,
          class BaseShifterType>
std::vector<storage::count_t> 
PdBG<BaseStorageType, BaseShifterType>
::insert_and_query_sequence(const std::string& sequence) {

    hashing::KmerIterator<extender_type> iter(sequence, &partitioner);
    std::vector<storage::count_t> counts(sequence.length() - _K + 1);

    size_t pos = 0;
    while(!iter.done()) {
        auto h = iter.next();
        counts[pos] = insert_and_query(h);
        ++pos;
    }

    return counts;
}


template <class BaseStorageType,
          class BaseShifterType>
std::vector<storage::count_t> 
PdBG<BaseStorageType, BaseShifterType>
::query_sequence(const std::string& sequence) {

    hashing::KmerIterator<extender_type> iter(sequence, &partitioner);
    std::vector<storage::count_t> counts(sequence.length() - _K + 1);

    size_t pos = 0;
    while(!iter.done()) {
        auto h = iter.next();
        counts[pos] = query(h);
        ++pos;
    }

    return counts;
}


template <class BaseStorageType,
          class BaseShifterType>
std::vector<storage::count_t> 
PdBG<BaseStorageType, BaseShifterType>
::query_sequence_rolling(const std::string& sequence) {

    hashing::KmerIterator<extender_type> iter(sequence, &partitioner);
    std::vector<storage::count_t> counts(sequence.length() - _K + 1);
    
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


template <class BaseStorageType,
          class BaseShifterType>
void 
PdBG<BaseStorageType, BaseShifterType>
::query_sequence(const std::string&             sequence,
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


template <class BaseStorageType,
          class BaseShifterType>
void 
PdBG<BaseStorageType, BaseShifterType>
::query_sequence(const std::string&             sequence,
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

template class PdBG<storage::BitStorage, hashing::FwdRollingShifter>;
template class PdBG<storage::BitStorage, hashing::CanRollingShifter>;
template class PdBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>;
template class PdBG<storage::SparseppSetStorage, hashing::CanRollingShifter>;
template class PdBG<storage::ByteStorage, hashing::FwdRollingShifter>;
template class PdBG<storage::ByteStorage, hashing::CanRollingShifter>;
template class PdBG<storage::NibbleStorage, hashing::FwdRollingShifter>;
template class PdBG<storage::NibbleStorage, hashing::CanRollingShifter>;
template class PdBG<storage::QFStorage, hashing::FwdRollingShifter>;
template class PdBG<storage::QFStorage, hashing::CanRollingShifter>;

}


