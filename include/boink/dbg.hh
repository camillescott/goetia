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
#include "storage.hh"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace boink {

using std::string;
using std::make_pair;

typedef std::pair<bool, bool> bit_pair_t;
typedef std::vector<bit_pair_t> bit_pair_vector_t;


template <class StorageType,
          class HashShifter>
class dBG : public KmerClient {
    StorageType S;
    HashShifter hasher;

public:

    const uint16_t N;
    const uint64_t max_table;
    const std::vector<uint64_t> sizes;

    typedef HashShifter shifter_type;
	typedef AssemblerMixin<dBG<StorageType, HashShifter>> assembler_type;
    typedef KmerIterator<HashShifter> kmer_iter_type;
    

    explicit dBG(uint16_t K, std::vector<uint64_t> storage_size) :
        KmerClient(K),
        S(storage_size),
        hasher(K),
        sizes(storage_size),
        N(storage_size.size()),
        max_table(*std::max_element(storage_size.begin(), storage_size.end())) {
        
    }

    explicit dBG(uint16_t K, uint64_t max_table, uint16_t N) :
        KmerClient(K),
        hasher(K),
        sizes(oxli::get_n_primes_near_x(N, max_table)),
        N(N),
        max_table(max_table),
        S(sizes) {

    }

    std::unique_ptr<dBG<StorageType, HashShifter>> clone() const {
        return std::make_unique<dBG<StorageType, HashShifter>>(_K, sizes);
    }

    hash_t hash(const string& kmer) const {
        return hasher.hash(kmer);
    }

    hash_t hash(const char * kmer) const {
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

    double estimated_fp() {
        double fp = n_occupied() / sizes[0];
        fp = pow(fp, N);
        return fp;
    }

    uint64_t add_sequence(const string& sequence,
                          std::vector<hash_t>& kmer_hashes,
                          std::vector<bool>& is_new) {
        KmerIterator<HashShifter> iter(sequence, _K);

        uint64_t n_consumed = 0;
        size_t pos = 0;
        bool kmer_consumed;
        while(!iter.done()) {
            hash_t h = iter.next();
            kmer_consumed = add(h);
            kmer_hashes.push_back(h);
            is_new.push_back(kmer_consumed);
            n_consumed += kmer_consumed;
            ++pos;
        }

        return n_consumed;
    }

    uint64_t add_sequence(const string& sequence,
                          std::set<hash_t>& new_kmers) {
        KmerIterator<HashShifter> iter(sequence, _K);

        uint64_t n_consumed = 0;
        size_t pos = 0;
        bool is_new;
        while(!iter.done()) {
            hash_t h = iter.next();
            is_new = add(h);
            if (is_new) {
                new_kmers.insert(h);
            }
            n_consumed += is_new;
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

    void get_counts(const string& sequence,
                        std::vector<count_t>& counts,
                        std::vector<hash_t>& hashes,
                        std::set<hash_t>& new_hashes) {

        KmerIterator<HashShifter> iter(sequence, _K);

        while(!iter.done()) {
            hash_t h = iter.next();
            count_t result = get(h);
            counts.push_back(result);
            hashes.push_back(h);
            if (result == 0) {
                new_hashes.insert(h);
            }
        }
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

    unique_ptr<KmerIterator<HashShifter>> get_hash_iter(const string& sequence) {
        return make_unique<KmerIterator<HashShifter>>(sequence, _K);
    }

    unique_ptr<assembler_type> get_assembler() {
        return make_unique<assembler_type>(this);
    }

};


typedef dBG<oxli::BitStorage, DefaultShifter> DefaultDBG;
}

#endif
