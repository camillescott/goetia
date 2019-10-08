/**
 * (c) Camille Scott, 2019
 * File   : hashshifter.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.07.2019
 */
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

#include <algorithm>
#include <cstring>
#include <deque>
#include <iterator>
#include <string>
#include <vector>

#include "boink/kmers/kmerclient.hh"
#include "boink/boink.hh"
#include "boink/hashing/alphabets.hh"
#include "boink/hashing/exceptions.hh"

#include "boink/ring_span.hpp"

namespace boink {
namespace hashing {


template <class Derived,
          class HashType = uint64_t>
class HashShifter : public kmers::KmerClient {
protected:

    char *                  kmer_buffer;
    nonstd::ring_span<char> kmer_window;

public:

    // the type to hash into. generally, uint64_t
    typedef HashType hash_type;

    // allowable alphabet
    const std::string& symbols;

    // Type for representing a neighbor hash with its prefix or suffix symbol
    struct shift_type {
        hash_type hash;
        char symbol;

        shift_type() : 
            hash(0),
            symbol('A') 
        {
        }

        shift_type(hash_type hash, char symbol) : 
            hash(hash),
            symbol(symbol)
        {
        }

        friend inline std::ostream&
        operator<<(std::ostream& os, const shift_type& shift) {
            os << "<shift_type symbol=" << shift.symbol << " hash=" << shift.hash << ">";
            return os;
        }
    };

    // A k-mer string and its hash value.
    struct kmer_type {
        hash_type hash;
        std::string kmer;

        kmer_type() :
            hash(0)
        {
        }

        kmer_type(const hash_type hash, const std::string kmer) :
            hash(hash),
            kmer(kmer)
        {
        }

        bool operator==(const kmer_type& other) const {
            return other.hash == this->hash;
        }

        friend inline std::ostream&
        operator<<(std::ostream& os, const kmer_type& kmer) {
            os << "<kmer_type kmer=" << kmer.kmer << " hash=" << kmer.hash << ">";
            return os;
        }

    };


    hash_type set_cursor(const std::string& sequence);

    hash_type set_cursor(const char * sequence);

    template <typename Iterator>
    hash_type set_cursor(const Iterator begin, const Iterator end);

    // shadowed by derived
    hash_type get() {
        return derived().get();
    }

    hash_type hash(const std::string& sequence) const {
        if (sequence.length() < _K) {
            throw SequenceLengthException("Sequence must at least length K");
        }
        _validate(sequence);
        return derived()._hash(sequence);
    }

    hash_type hash(const char * sequence) const {
        _validate(sequence);
        return derived()._hash(sequence);
    }

    // shadowed by derived impl
    std::vector<shift_type> gather_left() {
        return derived().gather_left();
    }

    hash_type shift_left(const char c) {
        _validate(c);
        hash_type h = derived().update_left(c);
        kmer_window.push_front(c);
        return h;
    }

    // shadowed by derived impl
    std::vector<shift_type> gather_right() {
        return derived().gather_right();
    }

    hash_type shift_right(const char c) {
        _validate(c);
        hash_type h = derived().update_right(c);
        kmer_window.push_back(c);
        return h;
    }

    const std::string get_cursor() const {
        return std::string(kmer_window.begin(),
                           kmer_window.end());
    }

    void get_cursor(std::deque<char>& d) const {

        for (auto symbol : kmer_window) {
            d.push_back(symbol);
        }
    }

    void get_cursor(kmer_type& result) {
        result.hash = get();
        result.kmer = get_cursor();
    }

private:

    HashShifter(const    std::string& start,
                uint16_t K,
                const    std::string& alphabet=DNA_SIMPLE)
        : KmerClient(K),
          kmer_buffer(new char[K]),
          kmer_window(kmer_buffer, kmer_buffer + K, kmer_buffer, K),
          initialized(false),
          symbols(alphabet)
    {
        load(start);
    }

    HashShifter(uint16_t K,
                const    std::string& alphabet=DNA_SIMPLE)
        : KmerClient(K),
          kmer_buffer(new char[K]),
          kmer_window(kmer_buffer, kmer_buffer + K, kmer_buffer, K),
          initialized(false),
          symbols(alphabet)
    {
    }

    ~HashShifter() {
        delete [] kmer_buffer;
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
        if(!this->is_valid(c)) {
            std::string msg("Invalid symbol: ");
            msg += c;
            throw InvalidCharacterException(msg.c_str());
        }
    }

    void _validate(const char * sequence) const {
        if (!this->is_valid(sequence)) {
            std::string msg("Invalid symbol in: ");
            msg += sequence;
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

    void load(const std::string& sequence) {
        //load(sequence.cbegin());
        for (uint16_t i = 0; i < this->_K; ++i)  {
            kmer_window.push_back(sequence[i]);
        }
    }

    void load(const char * sequence) {
        for (uint16_t i = 0; i < this->_K; ++i)  {
            kmer_window.push_back(sequence[i]);
        }
    }

    template <typename Iterator>
    void load(Iterator begin, Iterator end) {
        std::copy(begin, end, kmer_window);
    }

    template <typename Iterator>
    void load(Iterator begin) {
        for (uint16_t i = 0; i < this->_K; ++i)  {
            kmer_window.push_back(*begin);
            begin++;
        }
    }

};

} // hashing
} // boink

#endif
