/* hasher.hh -- k-mer hash functions
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef HASHING_HH
#define HASHING_HH

#include "boink.hh"
#include "oxli/alphabets.hh"
#include "oxli/kmer_hash.hh"

#include <deque>

namespace boink {

typedef uint64_t hash_t;
typedef std::pair<hash_t, hash_t> full_hash_t;


class KmerClient {
protected:
    const uint16_t _K;
public:
    explicit KmerClient(uint16_t K) : _K(K) {}
    const uint16_t K() const { return _K; }
};


template <class Derived,
          const std::string& Alphabet>
class HashShifter : public KmerClient {
public:
    static const std::string symbols;
    std::deque<char> symbol_deque;
    using KmerClient::KmerClient;

    HashShifter(const std::string& start, uint16_t K) :
        KmerClient(K) {

        reset(start);
    }
    
    hash_t reset(const std::string& sequence) {
        if (sequence.length() < _K) {
            throw BoinkException("Sequence must be length K");
        }
        symbol_deque.clear();
        symbol_deque.insert(symbol_deque.begin(), sequence.begin(), sequence.begin()+_K);
        return derived().reset();
    }

    hash_t reset() {
        return derived().reset();
    }

    hash_t hash() {
        return derived().hash();
    }

    hash_t hash(std::string& sequence) {
        return derived().hash(sequence);
    }

    std::vector<hash_t> shift_left() {
        return derived().shift_left();
    }

    hash_t shift_left(const char c) {
        symbol_deque.push_front(c);
        hash_t h = derived().shift_left(c);
        symbol_deque.pop_back();
        return h;
    }

    std::vector<hash_t> shift_right() {
        return derived().shift_right();
    }

    hash_t shift_right(const char c) {
        symbol_deque.push_back(c);
        hash_t h = derived().shift_right(c);
        symbol_deque.pop_front();
        return h;
    }

    std::string get_cursor() {
        return std::string(symbol_deque.begin(), symbol_deque.end());
    }

private:

    Derived& derived() {
        return *static_cast<Derived*>(this);
    }

};

template <const std::string& Alphabet>
class RollingHashShifter : public HashShifter<RollingHashShifter<Alphabet>,
                                              Alphabet> {
    CyclicHash<hash_t> hasher;

public:

    typedef HashShifter<RollingHashShifter<Alphabet>, Alphabet> BaseShifter;
    using BaseShifter::HashShifter;
    using BaseShifter::reset;

    RollingHashShifter(const std::string& start, uint16_t K) :
        BaseShifter(start, K), hasher(K) {}

    RollingHashShifter(uint16_t K) :
        BaseShifter(K), hasher(K) {}

    hash_t reset() {
        // This pattern "resets" a stack object
        // by calling its destructor and reallocating
        (&hasher)->~CyclicHash<hash_t>();
        new (&hasher) CyclicHash<hash_t>(this->_K);

        for (auto c : this->symbol_deque) {
            hasher.eat(c);
        }
        return hasher.hashvalue;
    }

    hash_t hash() {
        return hasher.hashvalue;
    }

    hash_t hash(std::string& sequence) {
        return _hash_cyclic_forward(sequence, this->_K);
    }

    hash_t shift_left(const char c) {
        hasher.reverse_update(c, this->symbol_deque.back());
        return hasher.hashvalue;
    }

    std::vector<hash_t> shift_right() {
        std::vector<hash_t> hashes;
        const char front = this->symbol_deque.front();
        for (auto symbol : Alphabet) {
            hasher.update(front, symbol);
            hashes.push_back(hasher.hashvalue);
            hasher.reverse_update(front, symbol);
        }
        return hashes;
    }

    hash_t shift_right(const char c) {
        hasher.update(this->symbol_deque.front(), c);
        return hash();
    }

    std::vector<hash_t> shift_left() {
        std::vector<hash_t> hashes;
        const char back = this->symbol_deque.back();
        for (auto symbol : Alphabet) {
            hasher.reverse_update(symbol, back);
            hashes.push_back(hasher.hashvalue);
            hasher.update(symbol, back);
        }

        return hashes;
    }
};

typedef RollingHashShifter<oxli::alphabets::DNA_SIMPLE> DefaultShifter;


template <class Shifter>
class KmerIterator : public KmerClient {
    const std::string _seq;
    unsigned int index;
    unsigned int length;
    bool _initialized;
    Shifter shifter;

public:
    KmerIterator(const std::string seq, uint16_t K) :
        KmerClient(K), _seq(seq), 
        index(0), _initialized(false), shifter(K) {
        
    }

    hash_t first() {
        _initialized = true;

        shifter.reset(_seq);
        index += 1;
        return shifter.hash();
    }

    hash_t next() {
        if (!_initialized) {
            return first();
        }

        if (done()) {
            throw BoinkException("past end of iterator");
        }

        shifter.shift_right(_seq[index + _K - 1]);
        index += 1;

        return shifter.hash();
    }

    bool done() const {
        return (index + _K > _seq.length());
    }

    unsigned int get_start_pos() const {
        if (!_initialized) { return 0; }
        return index - 1;
    }

    unsigned int get_end_pos() const {
        if (!_initialized) { return _K; }
        return index + _K - 1;
    }
};


/*
class FullRollingHasher {
    const char * _seq;
    const std::string _rev;
    const char _ksize;
    unsigned int index;
    unsigned int length;
    bool _initialized;
    oxli::CyclicHash<uint64_t> fwd_hasher;
    oxli::CyclicHash<uint64_t> bwd_hasher;

public:
    FullRollingHasher(const char * seq, unsigned char k) :
        _seq(seq), _rev(oxli::_revcomp(seq)), _ksize(k), index(0),
        _initialized(false), fwd_hasher(k), bwd_hasher(k)
    {
        length = strlen(_seq);
    };

    full_hash_t first() {
        _initialized = true;

        for (char i = 0; i < _ksize; ++i) {
            fwd_hasher.eat(*(_seq + i));
            bwd_hasher.eat(_rev[length - _ksize + i]);
        }
        index += 1;
        return std::make_pair(fwd_hasher.hashvalue, bwd_hasher.hashvalue);
    }

    full_hash_t next() {
        if (!_initialized) {
            return first();
        }

        if (done()) {
            throw oxli_exception("past end of iterator");
        }
        fwd_hasher.update(*(_seq + index - 1), *(_seq + index + _ksize - 1));

        // first argument is added, second is removed from the hash
        bwd_hasher.reverse_update(
          _rev[length - _ksize - index], _rev[length - index]);

        index += 1;

        return std::make_pair(fwd_hasher.hashvalue, bwd_hasher.hashvalue);
    }

    bool done() const {
        return (index + _ksize > length);
    }

    unsigned int get_start_pos() const {
        if (!_initialized) { return 0; }
        return index - 1;
    }

    unsigned int get_end_pos() const {
        if (!_initialized) { return _ksize; }
        return index + _ksize - 1;
    }
};
*/

};

#endif
