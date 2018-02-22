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

#include "hasher.hh"
#include "oxli/storage.hh"

#include <string>
#include <vector>

namespace boink {

typedef std::pair<bool, bool> bit_pair_t;
typedef std::vector<bit_pair_t> bit_pair_vector_t;

typedef uint8_t count_t;
typedef std::pair<uint8_t, uint8_t> full_count_t;


template <class StorageType,
          class HashShifter>
class dBG : public KmerClient {
    StorageType S;
    HashShifter hasher;

public:

    typedef HashShifter hasher_type;
    
    explicit dBG(uint16_t K, std::vector<uint64_t> storage_size) :
        KmerClient(K), S(storage_size), hasher(K) {

    }

    hash_t hash(const string& kmer) const {
        return hasher.hash(kmer);
    }

    /*
    full_hash_t hash_full(const string& kmer) const {
        hash_t fw, rc;
        oxli::_hash_cyclic(kmer, _K, fw, rc);
        return std::make_pair(fw, rc);
    }
    */

    bool add(const string& kmer) {
        return S.add(hash(kmer));
    }

    bool add(hash_t kmer) {
        return S.add(kmer);
    }

    /*
    bit_pair_t add_full(full_hash_t hashed) {
        return std::make_pair(S.add(hashed.first),
                              S.add(hashed.second));
    }

    bit_pair_t add_full(const string& kmer) {
        return add_full(hash_full(kmer));
    }
    */

    count_t get(const string& kmer) const {
        return S.get_count(hash(kmer));
    }

    count_t get(hash_t& hashed_kmer) const {
        return S.get_count(hashed_kmer);
    }

    /*
    full_count_t get(full_hash_t fully_hashed_kmer) const {
        return std::make_pair(get(fully_hashed_kmer.first),
                              get(fully_hashed_kmer.second));
    }

    full_count_t get_full(const string& kmer) const {
        return get(hash_full(kmer));
    }
    */

    std::vector<bool> add_sequence(const string& sequence);

    void save(string filename) const {
        S.save(filename, _K);
    }

    void load(string filename) {
        S.load(filename, _K);
    }
};


typedef dBG<oxli::BitStorage, DefaultShifter> DefaultDBG;

};

#endif
