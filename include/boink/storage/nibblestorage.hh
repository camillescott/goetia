/* nibblestorage.hh -- boink-modified oxli storage
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

#ifndef BOINK_NIBBLESTORAGE_HH
#define BOINK_NIBBLESTORAGE_HH

#include <cassert>
#include <cstring>
#include <array>
#include <memory>
#include <mutex>
#include <unordered_map>

#include "boink/hashing/hashing_types.hh"
#include "boink/storage/storage.hh"


namespace boink {
namespace storage {


/*
 * \class NibbleStorage
 *
 * \brief A A CountMin sketch implementation using 4bit counters.
 *
 * NibbleStorage is used to track counts of k-mers by Hashtable
 * and derived classes.  It contains 'n_tables' different tables of
 * 'tablesizes' entries. It allocates half a byte per table entry.
 *
 * Like other Storage classes, NibbleStorage manages setting the bits and
 * tracking statistics, as well as save/load, and not much else.
 *
 */
class NibbleStorage : public Storage
{
protected:
    // table size is measured in number of entries in the table, not in bytes
    std::vector<uint64_t> _tablesizes;
    size_t _n_tables;
    uint64_t _occupied_bins;
    uint64_t _n_unique_kmers;
    std::array<std::mutex, 32> mutexes;
    static constexpr uint8_t _max_count{15};
    byte_t ** _counts;

    // Compute index into the table, this retrieves the correct byte
    // which you then need to select the correct nibble from
    uint64_t _table_index(const hashing::hash_t k, const uint64_t tablesize) const
    {
        return (k % tablesize) / 2;
    }
    // Compute which half of the byte to use for this hash value
    uint8_t _mask(const hashing::hash_t k, const uint64_t tablesize) const
    {
        return (k%tablesize)%2 ? 15 : 240;
    }
    // Compute which half of the byte to use for this hash value
    uint8_t _shift(const hashing::hash_t k, const uint64_t tablesize) const
    {
        return (k%tablesize)%2 ? 0 : 4;
    }

public:
    NibbleStorage(uint64_t max_table, uint16_t N)
        : NibbleStorage(get_n_primes_near_x(N, max_table))
    {
    }

    NibbleStorage(const std::vector<uint64_t>& tablesizes) :
        _tablesizes{tablesizes},
        _n_tables(_tablesizes.size()),
        _occupied_bins{0},
        _n_unique_kmers{0}
    {
        // to allow more than 32 tables increase the size of mutex pool
        assert(_n_tables <= 32);
        _allocate_counters();
    }

    ~NibbleStorage()
    {
        if (_counts) {
            for (size_t i = 0; i < _n_tables; i++) {
                delete[] _counts[i];
                _counts[i] = NULL;
            }
            delete[] _counts;
            _counts = NULL;
            _n_tables = 0;
        }
    }

    std::unique_ptr<NibbleStorage> clone() const {
        return std::make_unique<NibbleStorage>(this->_tablesizes);
    }

    void _allocate_counters()
    {
        _counts = new byte_t*[_n_tables];

        for (size_t i = 0; i < _n_tables; i++) {
            const uint64_t tablesize = _tablesizes[i];
            const uint64_t tablebytes = tablesize / 2 + 1;

            _counts[i] = new byte_t[tablebytes];
            memset(_counts[i], 0, tablebytes);
        }
    }
    
    void reset()
    {
        for (unsigned int table_num = 0; table_num < _n_tables; table_num++) {
            uint64_t tablesize = _tablesizes[table_num];
            uint64_t tablebytes = tablesize / 2 + 1;
            memset(_counts[table_num], 0, tablebytes);
        }
    }

    inline const bool insert(hashing::hash_t khash)
    {
        bool is_new_kmer = false;

        for (unsigned int i = 0; i < _n_tables; i++) {
            MuxGuard g(mutexes[i]);
            byte_t* const table(_counts[i]);
            const uint64_t idx = _table_index(khash, _tablesizes[i]);
            const uint8_t mask = _mask(khash, _tablesizes[i]);
            const uint8_t shift = _shift(khash, _tablesizes[i]);
            const uint8_t current_count = (table[idx] & mask) >> shift;

            if (!is_new_kmer) {
                if (current_count == 0) {
                    is_new_kmer = true;

                    // track occupied bins in the first table only, as proxy
                    // for all.
                    if (i == 0) {
                        __sync_add_and_fetch(&_occupied_bins, 1);
                    }
                }
            }
            // if we have reached the maximum count stop incrementing the
            // counter. This avoids overflowing it.
            if (current_count == _max_count) {
                continue;
            }

            // increase count, no checking for overflow
            const uint8_t new_count = (current_count + 1) << shift;
            table[idx] = (table[idx] & ~mask) | (new_count & mask);
        }

        if (is_new_kmer) {
            __sync_add_and_fetch(&_n_unique_kmers, 1);
        }

        return is_new_kmer;
    }

    inline const count_t insert_and_query(hashing::hash_t khash)
    {
        if (insert(khash)) {
            return 1;
        }
        return query(khash);
    }

    // get the count for the given k-mer hash.
    inline const count_t query(hashing::hash_t khash) const
    {
        uint8_t min_count = _max_count; // bound count by maximum

        // get the minimum count across all tables
        for (unsigned int i = 0; i < _n_tables; i++) {
            const byte_t* table(_counts[i]);
            const uint64_t idx = _table_index(khash, _tablesizes[i]);
            const uint8_t mask = _mask(khash, _tablesizes[i]);
            const uint8_t shift = _shift(khash, _tablesizes[i]);
            const uint8_t the_count = (table[idx] & mask) >> shift;

            if (the_count < min_count) {
                min_count = the_count;
            }
        }
        return min_count;
    }

    // Accessors for protected/private table info members
    std::vector<uint64_t> get_tablesizes() const
    {
        return _tablesizes;
    }
    const size_t n_tables() const
    {
        return _n_tables;
    }
    const uint64_t n_unique_kmers() const
    {
        return _n_unique_kmers;
    }
    const uint64_t n_occupied() const
    {
        return _occupied_bins;
    }
    double estimated_fp() {
        double fp = n_occupied() / _tablesizes[0];
        fp = pow(fp, n_tables());
        return fp;
    }
    void save(std::string outfilename, uint16_t ksize);
    void load(std::string infilename, uint16_t& ksize);

    byte_t ** get_raw_tables()
    {
        return _counts;
    }
};


template<> 
struct is_probabilistic<NibbleStorage> { 
      static const bool value = true;
};

}
}

#endif
