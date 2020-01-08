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
    std::string& symbols;

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

        inline const hash_type value() const {
            return hash;
        }

        friend inline std::ostream&
        operator<<(std::ostream& os, const shift_type& shift) {
            os << "<shift_type symbol=" << shift.symbol << " hash=" << shift.value() << ">";
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

        inline const hash_type value() const {
            return hash;
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

    bool is_initialized() const {
        return initialized;
    }

private:

    HashShifter(const std::string& start,
                uint16_t           K,
                std::string&       alphabet=DNA_SIMPLE)
        : KmerClient(K),
          kmer_buffer(new char[K]),
          kmer_window(kmer_buffer, kmer_buffer + K, kmer_buffer, K),
          initialized(false),
          symbols(alphabet)
    {
        load(start);
    }

    HashShifter(uint16_t     K,
                std::string& alphabet=DNA_SIMPLE)
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
            std::string msg("HashShifter: Invalid symbol: ");
            msg += c;
            throw InvalidCharacterException(msg.c_str());
        }
    }

    void _validate(const char * sequence) const {
        if (!this->is_valid(sequence)) {
            std::string msg("HashShifter: Invalid symbol in: ");
            msg += sequence;
            throw InvalidCharacterException(msg.c_str());
        }
    }

    void _validate(const std::string& sequence) const {
        if (!is_valid(sequence)) {
            std::string msg("HashShifter: Invalid symbol in ");
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


template <class ShifterType>
class BiDirectionalShifter : public kmers::KmerClient {
protected:

    ShifterType fw_shifter, rc_shifter;

public:

    typedef ShifterType                      shifter_type;
    typedef typename ShifterType::hash_type  hash_type;

    struct shift_type {
        hash_type fw_hash, rc_hash;
        char fw_symbol, rc_symbol;
        bool sign;

        shift_type() :
            fw_hash(0),
            rc_hash(0),
            fw_symbol('A'),
            rc_symbol('A'),
            sign(true)
        {
        }

        shift_type(const hash_type fw_hash,
                   const hash_type rc_hash,
                   char fw_symbol,
                   char rc_symbol)
            : fw_hash(fw_hash),
              rc_hash(rc_hash),
              fw_symbol(fw_symbol),
              rc_symbol(rc_symbol),
              sign(fw_hash < rc_hash)
        {
        }

        inline const hash_type value() const {
            return sign ? fw_hash : rc_hash;
        }

    };

    struct kmer_type {
        hash_type   fw_hash, rc_hash;
        std::string kmer;
        bool        sign; // true if the string kmer is
                          // the canonical kmer

        kmer_type() :
            fw_hash(0),
            rc_hash(0),
            sign(true)
        {
        }

        kmer_type(const hash_type fw_hash,
                  const hash_type rc_hash,
                  const std::string kmer) :
            fw_hash(fw_hash),
            rc_hash(rc_hash),
            sign(fw_hash < rc_hash),
            kmer(kmer)
        {
        }

        const std::string rc() const {
            return revcomp(kmer);
        }

        inline const hash_type value() const {
            return sign ? fw_hash : rc_hash;
        }

        bool operator==(const kmer_type& other) const {
            return (this->sign == other.sign &&
                    this->fw_hash == other.fw_hash)
                   ||
                   (this->sign != other.sign &&
                    this->fw_hash == other.rc_hash);
        }

        operator hash_type() const {
            return sign ? fw_hash : rc_hash;
        }

        friend inline std::ostream&
        operator<<(std::ostream& os, const kmer_type& kmer) {
            os << "<kmer_type"
               << " fw=" << kmer.kmer
               << " rc=" << kmer.rc()
               << " sign=" << kmer.sign
               << " fwh=" << kmer.fw_hash
               << " rch=" << kmer.rc_hash
               << ">";
            return os;
        }
    };

    template<typename... Args>
    explicit BiDirectionalShifter(uint16_t K, Args&&... args)
        : fw_shifter(K, args...),
          rc_shifter(K, args...)
    {
        // okay this is hacky and inelegant but i don't feel
        // like major refactoring currently...
        rc_shifter.symbols = DNA_SIMPLE_CMP;
    }
 
    hash_type set_cursor(const std::string& sequence) {
        hash_type fw = fw_shifter.set_cursor(sequence);
        hash_type rc = rc_shifter.set_cursor(revcomp(sequence));
        return fw < rc ? fw : rc;
    }

    hash_type set_cursor(const char * sequence) {
        hash_type fw = fw_shifter.set_cursor(sequence);
        hash_type rc = rc_shifter.set_cursor(revcomp(fw_shifter.get_cursor()));
        return fw < rc ? fw : rc;
    }

    template <typename Iterator>
    hash_type set_cursor(const Iterator begin, const Iterator end) {
        hash_type fw = fw_shifter.set_cursor(begin, end);
        hash_type rc = rc_shifter.set_cursor(revcomp(fw_shifter.get_cursor()));
        return fw < rc ? fw : rc;
    }

    hash_type get() {
        hash_type fw = fw_shifter.get();
        hash_type rc = rc_shifter.get();
        return fw < rc ? fw : rc;
    }

    hash_type hash(const std::string& sequence) const {
        auto rc_kmer = revcomp(sequence.substr(0, this->K()));
        hash_type fw = fw_shifter._hash(sequence);
        hash_type rc = rc_shifter._hash(rc_kmer);
        //std::cout << fw << "  " << rc << std::endl;
        //std::cout << rc_kmer << std::endl;
        return fw < rc ? fw : rc;
    }

    std::vector<shift_type> gather_left() {
        auto fw_exts = fw_shifter.gather_left();
        auto rc_exts = rc_shifter.gather_right();

        std::vector<shift_type> result;
        for (size_t i = 0; i < fw_exts.size(); ++i) {
            auto& fw = fw_exts[i];
            auto& rc = rc_exts[i];
            result.emplace_back(fw.hash,
                                rc.hash,
                                fw.symbol,
                                rc.symbol);
        }
        
        return std::move(result);
    }

    hash_type shift_left(const char c);

    std::vector<shift_type> gather_right() {
        auto fw_exts = fw_shifter.gather_right();
        auto rc_exts = rc_shifter.gather_left();

        std::vector<shift_type> result;
        for (size_t i = 0; i < fw_exts.size(); ++i) {
            auto& fw = fw_exts[i];
            auto& rc = rc_exts[i];
            result.emplace_back(fw.hash,
                                rc.hash,
                                fw.symbol,
                                rc.symbol);
        }
        
        return std::move(result);
    }

    hash_type shift_right(const char c);

    const std::string get_cursor() const {
        return fw_shifter.get_cursor();
    }

    void get_cursor(std::deque<char>& d) const {
        fw_shifter.get_cursor(d);
    }

};

} // hashing
} // boink

#endif
