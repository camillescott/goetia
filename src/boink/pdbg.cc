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
    partitioner.reset_unikmers();

    return S->insert(partitioner.get(), partitioner.get_front_partition());
}


template <class BaseStorageType>
const bool 
PdBG<BaseStorageType>::insert(const hashing::PartitionedHash ph) {
    return S->insert(ph.first, ph.second);
}


template <class BaseStorageType>
const storage::count_t 
PdBG<BaseStorageType>::insert_and_query(const std::string& kmer) {
    partitioner.set_cursor(kmer);
    partitioner.reset_unikmers();

    return S->insert_and_query(partitioner.get(), partitioner.get_front_partition());
}


template <class BaseStorageType>
const storage::count_t 
PdBG<BaseStorageType>::insert_and_query(const hashing::PartitionedHash ph) {
    return S->insert_and_query(ph.first, ph.second);
}


template <class BaseStorageType>
const storage::count_t 
PdBG<BaseStorageType>::query(const std::string& kmer) {
    partitioner.set_cursor(kmer);
    partitioner.reset_unikmers();

    return S->query(partitioner.get(), partitioner.get_front_partition());
}


template <class BaseStorageType>
const storage::count_t 
PdBG<BaseStorageType>::query(const hashing::PartitionedHash ph) {
    return S->query(ph.first, ph.second);
}


template <class BaseStorageType>
uint64_t 
PdBG<BaseStorageType>::insert_sequence(const std::string&          sequence,
                         std::vector<hashing::hash_t>&  kmer_hashes,
                         std::vector<storage::count_t>& counts) {

    hashing::KmerIterator<UKHShifter> iter(sequence, &partitioner);

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


template <class BaseStorageType>
uint64_t 
PdBG<BaseStorageType>::insert_sequence(const std::string&      sequence,
                      std::set<hashing::hash_t>& new_kmers) {

    hashing::KmerIterator<UKHShifter> iter(sequence, &partitioner);

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


template <class BaseStorageType>
uint64_t 
PdBG<BaseStorageType>::insert_sequence(const std::string& sequence) {
    hashing::KmerIterator<UKHShifter> iter(sequence, &partitioner);

    uint64_t n_consumed = 0;
    while(!iter.done()) {
        hashing::PartitionedHash h = iter.next();
        n_consumed += insert(h);
    }

    return n_consumed;
}


template <class BaseStorageType>
uint64_t 
PdBG<BaseStorageType>::insert_sequence_rolling(const std::string& sequence) {
    hashing::KmerIterator<UKHShifter> iter(sequence, &partitioner);

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


template <class BaseStorageType>
std::vector<storage::count_t> 
PdBG<BaseStorageType>::insert_and_query_sequence(const std::string& sequence) {

    hashing::KmerIterator<UKHShifter> iter(sequence, &partitioner);
    std::vector<storage::count_t> counts(sequence.length() - _K + 1);

    size_t pos = 0;
    while(!iter.done()) {
        hashing::PartitionedHash h = iter.next();
        counts[pos] = insert_and_query(h);
        ++pos;
    }

    return counts;
}


template <class BaseStorageType>
std::vector<storage::count_t> 
PdBG<BaseStorageType>::query_sequence(const std::string& sequence) {

    hashing::KmerIterator<UKHShifter> iter(sequence, &partitioner);
    std::vector<storage::count_t> counts(sequence.length() - _K + 1);

    size_t pos = 0;
    while(!iter.done()) {
        hashing::PartitionedHash h = iter.next();
        counts[pos] = query(h);
        ++pos;
    }

    return counts;
}


template <class BaseStorageType>
std::vector<storage::count_t> 
PdBG<BaseStorageType>::query_sequence_rolling(const std::string& sequence) {

    hashing::KmerIterator<UKHShifter> iter(sequence, &partitioner);
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


template <class BaseStorageType>
void 
PdBG<BaseStorageType>::query_sequence(const std::string&             sequence,
                    std::vector<storage::count_t>& counts,
                    std::vector<hashing::hash_t>&  hashes) {

    hashing::KmerIterator<UKHShifter> iter(sequence, &partitioner);

    while(!iter.done()) {
        hashing::PartitionedHash h = iter.next();
        storage::count_t result = query(h);
        counts.push_back(result);
        hashes.push_back(h.first);
    }
}


template <class BaseStorageType>
void 
PdBG<BaseStorageType>::query_sequence(const std::string& sequence,
                    std::vector<storage::count_t>& counts,
                    std::vector<hashing::hash_t>& hashes,
                    std::set<hashing::hash_t>& new_hashes) {

    hashing::KmerIterator<UKHShifter> iter(sequence, &partitioner);

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

}

template class boink::PdBG<boink::storage::BitStorage>;
template class boink::PdBG<boink::storage::ByteStorage>;
template class boink::PdBG<boink::storage::NibbleStorage>;
template class boink::PdBG<boink::storage::QFStorage>;
template class boink::PdBG<boink::storage::SparseppSetStorage>;

