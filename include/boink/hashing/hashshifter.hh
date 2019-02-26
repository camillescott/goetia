/* hashshifter.hh -- base CRTP HashShifter class
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_HASHSHIFTER_HH
#define BOINK_HASHSHIFTER_HH

#include <deque>
#include <string>
#include <vector>

#include "boink/boink.hh"
#include "boink/hashing/alphabets.hh"
#include "boink/hashing/hashing_types.hh"

namespace boink {
namespace hashing {


class KmerClient {
protected:
    const uint16_t _K;
    explicit KmerClient(uint16_t K) : _K(K) {}

public:
    const uint16_t K() const { return _K; }
};


template <class Derived,
          const std::string& Alphabet = DNA_SIMPLE>
class HashShifter : public KmerClient {
public:

    typedef hashing::hash_t hash_type;

    static const std::string symbols;
    std::deque<char> symbol_deque;
    
    hash_t set_cursor(const std::string& sequence) {
        if (sequence.length() < _K) {
            throw SequenceLengthException("Sequence must at least length K");
        }
        if (!initialized) {
            load(sequence);
            derived().init();
        } else {
            for (auto cit = sequence.begin(); cit < sequence.begin() + _K; ++cit) {
                shift_right(*cit);
            }
        }
        return get();
    }

    hash_t set_cursor(const char * sequence) {
        // less safe! does not check length
        if(!initialized) {
            load(sequence);
            derived().init();
        } else {
            for (uint16_t i = 0; i < this->_K; ++i) {
                shift_right(sequence[i]);
            }
        }
        return get();
    }

    template <typename Iterator>
    hash_t set_cursor(const Iterator begin, const Iterator end) {
        Iterator _begin = begin, _end = end;
        if(!initialized) {
            load(begin, end);
            derived().init();
        } else {
            uint16_t l = 0;
            while (_begin != _end) {
                shift_right(*_begin);
                ++_begin;
                ++l;
            }
            if (l < this->_K) {
                throw SequenceLengthException("Sequence must at least length K");
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

    bool is_valid(const char * sequence) const {
        for (uint16_t i = 0; i < this->_K; ++i) {
            if(!is_valid(sequence[i])) {
                return false;
            }
        }
        return true;
    }

    void _validate(const char c) const {
        if (!this->is_valid(c)) {
            std::string msg("Invalid symbol: ");
            msg += c;
            throw InvalidCharacterException(msg.c_str());
        }
    }

    void _validate(const std::string& sequence) const {
        if (!is_valid(sequence)) {
            std::string msg("Invalid symbol in ");
            msg += sequence;
            msg += ", alphabet=";
            msg += symbols;
            throw InvalidCharacterException(msg.c_str());
        }
    }

    // shadowed by derived
    hash_t get() {
        return derived().get();
    }

    hash_t hash(const std::string& sequence) const {
        if (sequence.length() < _K) {
            throw SequenceLengthException("Sequence must at least length K");
        }
        _validate(sequence);
        return derived()._hash(sequence);
    }

    hash_t hash(const char * sequence) const {
        _validate(sequence);
        return derived()._hash(sequence);
    }

    // shadowed by derived impl
    std::vector<shift_t> gather_left() {
        return derived().gather_left();
    }

    hash_t shift_left(const char c) {
        _validate(c);
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
        _validate(c);
        symbol_deque.push_back(c);
        hash_t h = derived().update_right(c);
        symbol_deque.pop_front();
        return h;
    }

    std::string get_cursor() const {
        return std::string(symbol_deque.begin(), symbol_deque.end());
    }

    void get_cursor(std::deque<char>& d) const {
        d.insert(d.end(), symbol_deque.begin(), symbol_deque.end());
    }

private:

    HashShifter(const std::string& start, uint16_t K)
        : KmerClient(K), initialized(false)
    {
        load(start);
    }

    HashShifter(uint16_t K)
        : KmerClient(K), initialized(false)
    {
    }

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

    void load(const char * sequence) {
        symbol_deque.clear();
        for (uint16_t i = 0; i < this->_K; ++i) {
            symbol_deque.push_back(sequence[i]);
        }
    }

    template <typename Iterator>
    void load(const Iterator begin, const Iterator end) {
        symbol_deque.clear();
        symbol_deque.insert(symbol_deque.begin(),
                            begin,
                            end);
        if (symbol_deque.size() != this->_K) {
            throw SequenceLengthException("Sequence must be length K");
        }
    }

};


template<class Derived, const std::string& Alphabet>
const std::string HashShifter<Derived, Alphabet>::symbols = Alphabet;



} // hashing
} // boink

#endif
