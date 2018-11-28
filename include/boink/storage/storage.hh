/* storage.hh -- boink-modified oxli storage
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 *
 *** END BOINK LICENSE BLOCK
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

#ifndef BOINK_STORAGE_HH
#define BOINK_STORAGE_HH

#include <cmath>
#include <cassert>
#include <array>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <utility>
#include <vector>

#include "boink/hashing/hashing_types.hh"

namespace boink {
namespace storage {

using MuxGuard = std::lock_guard<std::mutex>;

typedef uint8_t                     byte_t;
typedef uint16_t                    count_t;
typedef std::pair<count_t, count_t> full_count_t;

typedef std::unordered_map<hashing::hash_t, count_t> KmerCountMap;

//
// base Storage class for hashtable-related storage of information in memory.
//

class Storage
{
protected:
    bool _supports_bigcount;
    bool _use_bigcount;

public:
    Storage() : _supports_bigcount(false), _use_bigcount(false) { } ;
    virtual ~Storage() { }
    virtual std::vector<uint64_t> get_tablesizes() const = 0;
    virtual const size_t n_tables() const = 0;
    virtual void save(std::string, uint16_t) = 0;
    virtual void load(std::string, uint16_t&) = 0;
    virtual const uint64_t n_occupied() const = 0;
    virtual const uint64_t n_unique_kmers() const = 0;
    virtual count_t test_and_set_bits( hashing::hash_t khash ) = 0;
    virtual bool add(hashing::hash_t khash) = 0;
    virtual const count_t get_count(hashing::hash_t khash) const = 0;
    virtual byte_t ** get_raw_tables() = 0;
    virtual void reset() = 0;

    void set_use_bigcount(bool b);
    bool get_use_bigcount();
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
