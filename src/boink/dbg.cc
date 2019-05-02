#include "boink/dbg.hh"

#include "boink/hashing/rollinghashshifter.hh"
#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/sparseppstorage.hh"

#include "boink/processors.hh"

#include <memory>


namespace boink {


template <class StorageType,
          class HashShifter>
uint64_t
dBG<StorageType, HashShifter>::insert_sequence(const std::string& sequence,
                                               std::vector<hashing::hash_t>&  kmer_hashes,
                                               std::vector<storage::count_t>& counts) {
    
    hashing::KmerIterator<HashShifter> iter(sequence, _K);

    uint64_t         n_consumed = 0;
    size_t           pos = 0;
    storage::count_t count;
    while(!iter.done()) {
        auto h = iter.next();
        count  = insert_and_query(h);

        kmer_hashes.push_back(h);
        counts.push_back(count);

        n_consumed += (count == 1);
        ++pos;
    }

    return n_consumed;
}


template <class StorageType,
          class HashShifter>
uint64_t
dBG<StorageType, HashShifter>::insert_sequence(const std::string&         sequence,
                                               std::set<hashing::hash_t>& new_kmers) {

    hashing::KmerIterator<HashShifter> iter(sequence, _K);

    uint64_t n_consumed = 0;
    size_t pos = 0;
    bool is_new;
    while(!iter.done()) {
        hashing::hash_t h = iter.next();
        if(insert(h)) {
            new_kmers.insert(h);
            n_consumed += is_new;
        }
        ++pos;
    }

    return n_consumed;
}


template <class StorageType,
          class HashShifter>
uint64_t
dBG<StorageType, HashShifter>::insert_sequence(const std::string& sequence) {
    
    hashing::KmerIterator<HashShifter> iter(sequence, _K);

    uint64_t n_consumed = 0;
    while(!iter.done()) {
        hashing::hash_t h = iter.next();
        n_consumed += insert(h);
    }

    return n_consumed;
}


template <class StorageType,
          class HashShifter>
std::vector<storage::count_t>
dBG<StorageType, HashShifter>::insert_and_query_sequence(const std::string& sequence) {

    hashing::KmerIterator<HashShifter> iter(sequence, _K);
    std::vector<storage::count_t> counts(sequence.length() - _K + 1);

    size_t pos = 0;
    while(!iter.done()) {
        hashing::hash_t h = iter.next();
        counts[pos] = S->insert_and_query(h);
        ++pos;
    }

    return counts;
}


template <class StorageType,
          class HashShifter>
std::vector<storage::count_t>
dBG<StorageType, HashShifter>::query_sequence(const std::string& sequence) {

    hashing::KmerIterator<HashShifter> iter(sequence, _K);
    std::vector<storage::count_t> counts(sequence.length() - _K + 1);

    size_t pos = 0;
    while(!iter.done()) {
        hashing::hash_t h = iter.next();
        counts[pos] = query(h);
        ++pos;
    }

    return counts;
}


template <class StorageType,
          class HashShifter>
void
dBG<StorageType, HashShifter>::query_sequence(const std::string&             sequence,
                                              std::vector<storage::count_t>& counts,
                                              std::vector<hashing::hash_t>&  hashes) {

    hashing::KmerIterator<HashShifter> iter(sequence, _K);

    while(!iter.done()) {
        hashing::hash_t h = iter.next();
        storage::count_t result = query(h);
        counts.push_back(result);
        hashes.push_back(h);
    }
}


template <class StorageType,
          class HashShifter>
void
dBG<StorageType, HashShifter>::query_sequence(const std::string& sequence,
                                              std::vector<storage::count_t>& counts,
                                              std::vector<hashing::hash_t>& hashes,
                                              std::set<hashing::hash_t>& new_hashes) {

	hashing::KmerIterator<HashShifter> iter(sequence, _K);

	while(!iter.done()) {
		hashing::hash_t h = iter.next();
		auto result = query(h);
		if (result == 0) {
			new_hashes.insert(h);
		}
		counts.push_back(result);
		hashes.push_back(h);
	}
}

}

template class boink::dBG<boink::storage::BitStorage,
                          boink::hashing::RollingHashShifter>;
template class boink::dBG<boink::storage::ByteStorage,
                          boink::hashing::RollingHashShifter>;
template class boink::dBG<boink::storage::NibbleStorage,
                          boink::hashing::RollingHashShifter>;
template class boink::dBG<boink::storage::QFStorage,
                          boink::hashing::RollingHashShifter>;
template class boink::dBG<boink::storage::SparseppSetStorage,
                          boink::hashing::RollingHashShifter>;


template class boink::InserterProcessor<boink::dBG<boink::storage::BitStorage,
                                                   boink::hashing::RollingHashShifter>>;
template class boink::InserterProcessor<boink::dBG<boink::storage::ByteStorage,
                                                   boink::hashing::RollingHashShifter>>;
template class boink::InserterProcessor<boink::dBG<boink::storage::NibbleStorage,
                                                   boink::hashing::RollingHashShifter>>;
template class boink::InserterProcessor<boink::dBG<boink::storage::QFStorage,
                                                   boink::hashing::RollingHashShifter>>;
template class boink::InserterProcessor<boink::dBG<boink::storage::SparseppSetStorage,
                                                   boink::hashing::RollingHashShifter>>;

