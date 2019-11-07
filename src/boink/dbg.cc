/**
 * (c) Camille Scott, 2019
 * File   : dbg.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */

#include "boink/dbg.hh"

#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/sparseppstorage.hh"

#include "boink/processors.hh"

#include <memory>


namespace boink {

/*
 * CONSTRUCTORS
 */ 

template<class StorageType,
         class ShifterType>
dBG<StorageType, ShifterType>

::dBG(std::shared_ptr<StorageType> S, ShifterType& hasher)
    : KmerClient(hasher.K()),
      hasher(hasher),
      S(S)
{
}

template<class StorageType,
         class ShifterType>
template<typename... Args>
dBG<StorageType, ShifterType>

::dBG(std::shared_ptr<StorageType> S, uint16_t K, Args&&... args)
    : KmerClient(K),
      hasher(K, std::forward<Args>(args)...),
      S(S)
{
}

/*
 * FACTORIES
 */

template<class StorageType,
         class ShifterType>
std::shared_ptr<dBG<StorageType, ShifterType>>
dBG<StorageType, ShifterType>

::build(std::shared_ptr<StorageType> S, ShifterType& hasher) {
    return std::make_shared<dBG<StorageType, ShifterType>>(S, hasher);
}


template<class StorageType,
         class ShifterType>
template<typename... Args>
std::shared_ptr<dBG<StorageType, ShifterType>>
dBG<StorageType, ShifterType>

::build(std::shared_ptr<StorageType> S, uint16_t K, Args&&... args) {
    return std::make_shared<dBG<StorageType, ShifterType>>(S, K, std::forward<Args>(args)...);
}


template<class StorageType,
         class ShifterType>
std::shared_ptr<dBG<StorageType, ShifterType>>
dBG<StorageType, ShifterType>

::clone() {
    return std::make_shared<dBG<StorageType, ShifterType>>(S->clone(), hasher);
}


/*
 * METHODS
 */

template<class StorageType,
         class ShifterType>
typename ShifterType::hash_type
dBG<StorageType, ShifterType>

::hash(const std::string& kmer) {
    hasher.set_cursor(kmer);
    return hasher.get();
}


template<class StorageType,
         class ShifterType>
typename ShifterType::hash_type
dBG<StorageType, ShifterType>

::hash(const char * kmer) {
    hasher.set_cursor(kmer);
    return hasher.get();
}


//
// query
//

template<class StorageType,
         class ShifterType>
const storage::count_t
dBG<StorageType, ShifterType>

::query(hash_type hashed_kmer) const {
    return S->query(hashed_kmer);
}


template<class StorageType,
         class ShifterType>
const storage::count_t
dBG<StorageType, ShifterType>

::query(const std::string& kmer) {
    return S->query(hash(kmer));
}

//
// neighbor searchings
//

template <class StorageType,
          class ShifterType>
std::vector<typename ShifterType::shift_type> 
dBG<StorageType, ShifterType>

::left_neighbors(const std::string& root) {
    typename traversal_type::dBG traverser(hasher);
    traverser.set_cursor(root);
    return traverser.filter_nodes(this, traverser.gather_left());
}


template <class StorageType,
          class ShifterType>
std::vector<typename ShifterType::shift_type>
dBG<StorageType, ShifterType>

::right_neighbors(const std::string& root) {
    typename traversal_type::dBG traverser(hasher);
    traverser.set_cursor(root);
    return traverser.filter_nodes(this, traverser.gather_right());
}


template <class StorageType,
          class ShifterType>
std::pair<std::vector<typename ShifterType::shift_type>,
          std::vector<typename ShifterType::shift_type>>
dBG<StorageType, ShifterType>

::neighbors(const std::string& root) {
    typename traversal_type::dBG traverser(hasher);
    traverser.set_cursor(root);
    auto lfiltered = traverser.filter_nodes(this, traverser.gather_left());
    auto rfiltered = traverser.filter_nodes(this, traverser.gather_right());
    return std::make_pair(lfiltered, rfiltered);
}

//
// inserts
//

template <class StorageType,
          class ShifterType>
uint64_t
dBG<StorageType, ShifterType>::insert_sequence(const std::string& sequence,
                                               std::vector<hash_type>&  kmer_hashes,
                                               std::vector<storage::count_t>& counts) {
    
    hashing::KmerIterator<ShifterType> iter(sequence, &hasher);

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
          class ShifterType>
uint64_t
dBG<StorageType, ShifterType>::insert_sequence(const std::string&         sequence,
                                               std::set<hash_type>& new_kmers) {

    hashing::KmerIterator<ShifterType> iter(sequence, &hasher);

    uint64_t n_consumed = 0;
    size_t pos = 0;
    bool is_new;
    while(!iter.done()) {
        hash_type h = iter.next();
        if(insert(h)) {
            new_kmers.insert(h);
            n_consumed += is_new;
        }
        ++pos;
    }

    return n_consumed;
}


template <class StorageType,
          class ShifterType>
uint64_t
dBG<StorageType, ShifterType>::insert_sequence(const std::string& sequence) {
    
    hashing::KmerIterator<ShifterType> iter(sequence, &hasher);

    uint64_t n_consumed = 0;
    while(!iter.done()) {
        hash_type h = iter.next();
        n_consumed += insert(h);
    }

    return n_consumed;
}


template <class StorageType,
          class ShifterType>
std::vector<storage::count_t>
dBG<StorageType, ShifterType>::insert_and_query_sequence(const std::string& sequence) {

    hashing::KmerIterator<ShifterType> iter(sequence, &hasher);
    std::vector<storage::count_t> counts(sequence.length() - _K + 1);

    size_t pos = 0;
    while(!iter.done()) {
        hash_type h = iter.next();
        counts[pos] = S->insert_and_query(h);
        ++pos;
    }

    return counts;
}


template <class StorageType,
          class ShifterType>
std::vector<storage::count_t>
dBG<StorageType, ShifterType>::query_sequence(const std::string& sequence) {

    hashing::KmerIterator<ShifterType> iter(sequence, &hasher);
    std::vector<storage::count_t> counts(sequence.length() - _K + 1);

    size_t pos = 0;
    while(!iter.done()) {
        hash_type h = iter.next();
        counts[pos] = query(h);
        ++pos;
    }

    return counts;
}


template <class StorageType,
          class ShifterType>
void
dBG<StorageType, ShifterType>::query_sequence(const std::string&             sequence,
                                              std::vector<storage::count_t>& counts,
                                              std::vector<hash_type>&  hashes) {

    hashing::KmerIterator<ShifterType> iter(sequence, &hasher);

    while(!iter.done()) {
        hash_type h = iter.next();
        storage::count_t result = query(h);
        counts.push_back(result);
        hashes.push_back(h);
    }
}


template <class StorageType,
          class ShifterType>
void
dBG<StorageType, ShifterType>::query_sequence(const std::string& sequence,
                                              std::vector<storage::count_t>& counts,
                                              std::vector<hash_type>& hashes,
                                              std::set<hash_type>& new_hashes) {

	hashing::KmerIterator<ShifterType> iter(sequence, &hasher);

	while(!iter.done()) {
		hash_type h = iter.next();
		auto result = query(h);
		if (result == 0) {
			new_hashes.insert(h);
		}
		counts.push_back(result);
		hashes.push_back(h);
	}
}


template<class StorageType,
         class ShifterType>
std::shared_ptr<hashing::KmerIterator<ShifterType>>
dBG<StorageType, ShifterType>::get_hash_iter(const std::string& sequence) {
    // Return a KmerIterator with a *copy* of the hasher
    return std::make_shared<hashing::KmerIterator<ShifterType>>(sequence, hasher);
}


template<class StorageType,
         class ShifterType>
ShifterType
dBG<StorageType, ShifterType>::get_hasher() {
    return ShifterType(hasher);
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


template class boink::dBG<boink::storage::BitStorage,
                          boink::hashing::UKHS::LazyShifter>;
template class boink::dBG<boink::storage::ByteStorage,
                          boink::hashing::UKHS::LazyShifter>;
template class boink::dBG<boink::storage::NibbleStorage,
                          boink::hashing::UKHS::LazyShifter>;
template class boink::dBG<boink::storage::QFStorage,
                          boink::hashing::UKHS::LazyShifter>;
template class boink::dBG<boink::storage::SparseppSetStorage,
                          boink::hashing::UKHS::LazyShifter>;


void boink::test_dbg() {
    //auto g = boink::dBG<boink::storage::SparseppSetStorage,
    //                    boink::hashing::RollingHashShifter>::build(std::make_tuple(5),
    //                                                               std::make_tuple(0));
}
