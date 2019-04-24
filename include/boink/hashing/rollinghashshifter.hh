/* rollinghashshifter.hh -- k-mer hash functions
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_ROLLINGHASHSHIFTER_HH
#define BOINK_ROLLINGHASHSHIFTER_HH

#include "boink/boink.hh"
#include "boink/hashing/alphabets.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/hashshifter.hh"

#include "boink/hashing/rollinghash/cyclichash.h"

namespace boink {
namespace hashing {


class RollingHashShifter : public HashShifter<RollingHashShifter> {
protected:
    typedef HashShifter<RollingHashShifter> BaseShifter;

    CyclicHash<hash_t> hasher;
    using BaseShifter::_K;

public:
    using BaseShifter::symbols;

    //using BaseShifter::HashShifter;
    typedef hash_t hash_type;

    RollingHashShifter(const std::string& start,
                       uint16_t K)
        : BaseShifter(start, K), hasher(K)
    {    
        init();
    }

    RollingHashShifter(uint16_t K)
        : BaseShifter(K),
          hasher(K)
    {
    }

    RollingHashShifter(const RollingHashShifter& other)
        : BaseShifter(other.K()),
          hasher(other.K())
    {
        this->load(other.get_cursor());
        init();
    }

    void init() {
        if (this->initialized) {
            return;
        }
        for (auto c : this->kmer_window) {
            this->_validate(c);
            hasher.eat(c);
        }
        this->initialized = true;
    }

    hash_t get() {
        return hasher.hashvalue;
    }

    hash_t _hash(const std::string& sequence) const {
        return hash_cyclic(sequence, this->_K);
    }

    hash_t _hash(const char * sequence) const {
        CyclicHash<hash_t> tmp_hasher(this->_K);
        for (uint16_t i = 0; i < this->_K; ++i) {
            tmp_hasher.eat(sequence[i]);
        }
        return tmp_hasher.hashvalue;
    }

    hash_t update_left(const char c) {
        hasher.reverse_update(c, this->kmer_window.back());
        return get();
    }

    hash_t update_right(const char c) {
        hasher.update(this->kmer_window.front(), c);
        return get();
    }

    std::vector<shift_t> gather_left() {
        std::vector<shift_t> hashes;
        const char back = this->kmer_window.back();
        for (auto symbol : symbols) {
            hasher.reverse_update(symbol, back);
            shift_t result(hasher.hashvalue, symbol);
            hashes.push_back(result);
            hasher.update(symbol, back);
        }

        return hashes;
    }

    std::vector<shift_t> gather_right() {
        std::vector<shift_t> hashes;
        const char front = this->kmer_window.front();
        for (auto symbol : symbols) {
            hasher.update(front, symbol);
            hashes.push_back(shift_t(hasher.hashvalue, symbol));
            hasher.reverse_update(front, symbol);
        }
        return hashes;
    }
};


} // hashing
} // boink

#endif
