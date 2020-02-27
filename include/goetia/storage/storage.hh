/**
 * (c) Camille Scott, 2019
 * File   : storage.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019

 *** END GOETIA LICENSE BLOCK
 *
 * This file is part of khmer, https://github.com/dib-lab/khmer/, and is
 * Copyright (C) 2016, The Regents of the University of California.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of the University of California nor the names
 *       of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written
 *       permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * LICENSE (END)
 * 
 * Contact: khmer-project@idyll.org
 */

#ifndef GOETIA_STORAGE_HH
#define GOETIA_STORAGE_HH

#include <cmath>
#include <cassert>
#include <array>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <utility>
#include <vector>

#include "goetia/goetia.hh"


#   define MAX_BIGCOUNT 65535
#   define SAVED_SIGNATURE "OXLI"
#   define SAVED_FORMAT_VERSION 4
#   define SAVED_COUNTING_HT 1
#   define SAVED_HASHBITS 2
#   define SAVED_TAGS 3
#   define SAVED_STOPTAGS 4
#   define SAVED_SUBSET 5
#   define SAVED_LABELSET 6
#   define SAVED_SMALLCOUNT 7
#   define SAVED_QFCOUNT 8


namespace goetia {
namespace storage {

using MuxGuard = std::lock_guard<std::mutex>;

typedef uint8_t                     byte_t;
typedef int16_t                     count_t;
typedef std::pair<count_t, count_t> full_count_t;

template<class Storage> 
struct is_probabilistic { 
    static const bool value = false;
};

template<class Storage>
struct is_counting {
    static const bool value = false;
};

//
// base Storage class for hashtable-related storage of information in memory.
//

template<typename ValueType>
class Storage
{
protected:
    bool _supports_bigcount;
    bool _use_bigcount;

public:

    typedef ValueType value_type;
    typedef std::unordered_map<value_type, count_t> CountMap;

    Storage() : _supports_bigcount(false), _use_bigcount(false) { } ;
    virtual ~Storage() { }

    //virtual std::vector<uint64_t> get_tablesizes() const = 0;
    //virtual const size_t n_tables() const = 0;

    virtual void save(std::string, uint16_t) = 0;
    virtual void load(std::string, uint16_t&) = 0;

    virtual const uint64_t n_occupied() const = 0;
    virtual const uint64_t n_unique_kmers() const = 0;

    virtual const bool    insert(value_type khash ) = 0;
    virtual const count_t insert_and_query(value_type khash) = 0;
    virtual const count_t query(value_type khash) const = 0;

    virtual byte_t ** get_raw_tables() = 0;
    virtual void reset() = 0;

    void set_use_bigcount(bool b)
    {
        if (!_supports_bigcount) {
            throw GoetiaException("bigcount is not supported for this storage.");
        }
        _use_bigcount = b;
    }

    bool get_use_bigcount()
    {
        return _use_bigcount;
    }

    virtual void serialize(std::ofstream& out) {}
};


inline bool is_prime(uint64_t n)
{
    if (n < 2) {
        return false;
    }
    if (n == 2) {
        return true;
    }
    if (n % 2 == 0) {
        return false;
    }
    for (unsigned long long i=3; i < sqrt(n) + 1; i += 2) {
        if (n % i == 0) {
            return false;
        }
    }
    return true;
}


inline std::vector<uint64_t> get_n_primes_near_x(uint32_t n, uint64_t x)
{
    std::vector<uint64_t> primes;
    if (x == 1) {
        primes.push_back(1);
        return primes;
    }

    uint64_t i = x - 1;
    if (i % 2 == 0) {
        i--;
    }
    while (primes.size() != n) {
        if (is_prime(i)) {
            primes.push_back(i);
        }
        if (i == 1) {
            break;
        }
        i -= 2;
    }

    // might return < n primes if x is too small
    return primes;
}


}
}

#endif 
