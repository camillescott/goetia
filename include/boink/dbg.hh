/* dbg.hh -- Generic de Bruijn graph
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BDBG_HH
#define BDBG_HH

#include "hashing.hh"
#include "assembly.hh"
#include "oxli/storage.hh"
#include "oxli/hashtable.hh"

#include <string>
#include <vector>

namespace boink {

using std::string;
using std::make_pair;
using std::shared_ptr;

typedef std::pair<bool, bool> bit_pair_t;
typedef std::vector<bit_pair_t> bit_pair_vector_t;

typedef uint8_t count_t;
typedef std::pair<uint8_t, uint8_t> full_count_t;


template <class StorageType,
          class HashShifter>
class dBG : public KmerClient,
		    public std::enable_shared_from_this<dBG<StorageType, HashShifter>> {
    StorageType S;
    HashShifter hasher;

public:

    typedef HashShifter shifter_type;
	typedef AssemblerMixin<dBG<StorageType, HashShifter>> assembler_type;
    

    explicit dBG(uint16_t K, std::vector<uint64_t> storage_size) :
        KmerClient(K),
        S(storage_size),
        hasher(K) {
        
    }

    explicit dBG(uint16_t K, uint64_t max_table, uint16_t N) :
        dBG(K, oxli::get_n_primes_near_x(N, max_table)) {

    }

	shared_ptr<dBG<StorageType, HashShifter>> get_ptr() {
		return this->shared_from_this();
	}

    hash_t hash(const string& kmer) const {
        return hasher.hash(kmer);
    }

    bool add(const string& kmer) {
        return S.add(hash(kmer));
    }

    bool add(hash_t kmer) {
        return S.add(kmer);
    }

    const count_t get(const string& kmer) const {
        return S.get_count(hash(kmer));
    }

    const count_t get(hash_t hashed_kmer) const {
        return S.get_count(hashed_kmer);
    }

    uint64_t n_unique() const {
        return S.n_unique_kmers();
    }

    uint64_t n_occupied() const {
        return S.n_occupied();
    }

    uint8_t ** get_raw() const {
        return S.get_raw_tables();
    }

    uint64_t add_sequence(const string& sequence,
                          std::vector<bool>& consumed) {
        KmerIterator<HashShifter> iter(sequence, _K);

        uint64_t n_consumed = 0;
        size_t pos = 0;
        bool kmer_consumed;
        while(!iter.done()) {
            hash_t h = iter.next();
            kmer_consumed = add(h);
            consumed.push_back(kmer_consumed);
            n_consumed += kmer_consumed;
            ++pos;
        }

        return n_consumed;
    }

    uint64_t add_sequence(const string& sequence) {
        KmerIterator<HashShifter> iter(sequence, _K);

        uint64_t n_consumed = 0;
        while(!iter.done()) {
            hash_t h = iter.next();
            n_consumed += add(h);
        }

        return n_consumed;
    }

    std::vector<count_t> get_counts(const string& sequence) {
        KmerIterator<HashShifter> iter(sequence, _K);
        std::vector<count_t> counts(sequence.length() - _K + 1);

        size_t pos = 0;
        while(!iter.done()) {
            hash_t h = iter.next();
            counts[pos] = get(h);
            ++pos;
        }

        return counts;
    }

    void save(string filename) {
        S.save(filename, _K);
    }

    void load(string filename) {
        unsigned char ksize = _K;
        S.load(filename, ksize);
    }

    void reset() {
        S.reset();
    }

    shared_ptr<KmerIterator<HashShifter>> get_hash_iter(const string& sequence) {
        return make_shared<KmerIterator<HashShifter>>(sequence, _K);
    }

    shared_ptr<assembler_type> get_assembler() {
        return make_shared<assembler_type>(this->get_ptr());
    }

};


typedef dBG<oxli::BitStorage, DefaultShifter> DefaultDBG;
}

#endif
