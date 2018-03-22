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

# ifdef DEBUG_HASHING
#   define pdebug(x) do { std::cerr << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif


namespace boink {

typedef uint64_t hash_t;
typedef std::pair<hash_t, hash_t> full_hash_t;

struct shift_t {
    hash_t hash;
    char symbol;

    shift_t() : 
        hash(0),
        symbol('A') {
    
    }

    shift_t(hash_t hash, char symbol) : 
        hash(hash),
        symbol(symbol) {
    
    }
};


class KmerClient {
protected:
    const uint16_t _K;
public:
    explicit KmerClient(uint16_t K) : _K(K) {}
    const uint16_t K() const { return _K; }
};


template <class Derived,
          const std::string& Alphabet = oxli::alphabets::DNA_SIMPLE>
class HashShifter : public KmerClient {
public:
    static const std::string symbols;
    std::deque<char> symbol_deque;
    
    hash_t set_cursor(const std::string& sequence) {
        if (sequence.length() < _K) {
            throw BoinkException("Sequence must be length K");
        }
        if (!initialized) {
            load(sequence);
            derived().init();
        } else {
            for (auto c : sequence) {
                shift_right(c);
            }
        }
        return get();
    }

    bool is_valid(const char c) const {
        for (auto symbol : symbols) {
            if (c == symbol) {
                return true;
            }
        }
        return false;
    }

    bool is_valid(const std::string& sequence) const {
        for (auto c : sequence) {
            if(!is_valid(c)) {
                return false;
            }
        }
        return true;
    }

    // shadowed by derived
    hash_t get() {
        return derived().get();
    }

    hash_t hash(const std::string& sequence) const {
        if (!is_valid(sequence)) {
            throw BoinkException("Invalid symbol.");
        }
        return derived()._hash(sequence);
    }

    // shadowed by derived impl
    std::vector<shift_t> gather_left() {
        return derived().gather_left();
    }

    hash_t shift_left(const char c) {
        if (!is_valid(c)) {
            throw BoinkException("Invalid symbol.");
        }
        symbol_deque.push_front(c);
        hash_t h = derived().update_left(c);
        symbol_deque.pop_back();
        return h;
    }

    // shadowed by derived impl
    std::vector<shift_t> gather_right() {
        return derived().gather_right();
    }

    hash_t shift_right(const char c) {
        if (!is_valid(c)) {
            throw BoinkException("Invalid symbol.");
        }
        symbol_deque.push_back(c);
        hash_t h = derived().update_right(c);
        symbol_deque.pop_front();
        return h;
    }

    std::string get_cursor() {
        return std::string(symbol_deque.begin(), symbol_deque.end());
    }

    void get_cursor(std::deque<char>& d) {
        d.insert(d.end(), symbol_deque.begin(), symbol_deque.end());
    }

private:

    HashShifter(const std::string& start, uint16_t K) :
        KmerClient(K), initialized(false) {

        load(start);
    }

    HashShifter(uint16_t K) :
        KmerClient(K), initialized(false) {}

    friend Derived;

    Derived& derived() {
        return *static_cast<Derived*>(this);
    }

    const Derived& derived() const {
        return *static_cast<const Derived*>(this);
    }

protected:

    bool initialized;

    void load(const std::string& sequence) {
        symbol_deque.clear();
        symbol_deque.insert(symbol_deque.begin(),
                            sequence.begin(),
                            sequence.begin()+_K); 
    }

};

template<class Derived, const std::string& Alphabet>
const std::string HashShifter<Derived, Alphabet>::symbols = Alphabet;



template <const std::string& Alphabet = oxli::alphabets::DNA_SIMPLE>
class RollingHashShifter : public HashShifter<RollingHashShifter<Alphabet>,
                                              Alphabet> {
    CyclicHash<hash_t> hasher;

public:

    typedef HashShifter<RollingHashShifter<Alphabet>, Alphabet> BaseShifter;
    //using BaseShifter::HashShifter;

    RollingHashShifter(const std::string& start, uint16_t K) :
        BaseShifter(start, K), hasher(K) {
        
        init();
    }

    RollingHashShifter(uint16_t K) :
        BaseShifter(K), hasher(K) {}

    void init() {
        if (this->initialized) {
            return;
        }
        for (auto c : this->symbol_deque) {
            if (!this->is_valid(c)) {
                throw BoinkException("Invalid symbol.");
            }
            hasher.eat(c);
        }
        this->initialized = true;
    }

    hash_t get() {
        return hasher.hashvalue;
    }

    hash_t _hash(const std::string& sequence) const {
        return oxli::_hash_cyclic_forward(sequence, this->_K);
    }

    hash_t update_left(const char c) {
        hasher.reverse_update(c, this->symbol_deque.back());
        return get();
    }

    std::vector<shift_t> gather_right() {
        std::vector<shift_t> hashes;
        const char front = this->symbol_deque.front();
        for (auto symbol : Alphabet) {
            hasher.update(front, symbol);
            hashes.push_back(shift_t(hasher.hashvalue, symbol));
            hasher.reverse_update(front, symbol);
        }
        return hashes;
    }

    hash_t update_right(const char c) {
        hasher.update(this->symbol_deque.front(), c);
        return get();
    }

    std::vector<shift_t> gather_left() {
        std::vector<shift_t> hashes;
        const char back = this->symbol_deque.back();
        for (auto symbol : Alphabet) {
            hasher.reverse_update(symbol, back);
            shift_t result(hasher.hashvalue, symbol);
            hashes.push_back(result);
            hasher.update(symbol, back);
        }

        return hashes;
    }
};

typedef RollingHashShifter<oxli::alphabets::DNA_SIMPLE> DefaultShifter;


template <class ShifterType>
class KmerIterator : public KmerClient {
    const std::string _seq;
    unsigned int index;
    unsigned int length;
    bool _initialized;
    ShifterType shifter;

public:
    KmerIterator(const std::string seq, uint16_t K) :
        KmerClient(K), _seq(seq), 
        index(0), _initialized(false), shifter(K) {

        if (_seq.length() < _K) {
            throw BoinkException("Sequence must have length >= K");
        }
        
    }

    hash_t first() {
        _initialized = true;

        shifter.set_cursor(_seq);
        index += 1;
        return shifter.get();
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

        return shifter.get();
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
