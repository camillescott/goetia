/* bitstorage.hh -- boink-modified oxli storage
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

#ifndef BOINK_BITSTORAGE_HH
#define BOINK_BITSTORAGE_HH

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
 * \class BitStorage
 *
 * \brief A Bloom filter implementation.
 *
 * BitStorage is used to track presence/absence of k-mers by Hashtable
 * and derived classes.  It contains 'n_tables' different tables of
 * bitsizes specified in 'tablesizes' (so 1/8 for bytesizes).
 *
 * Like other Storage classes, BitStorage manages setting the bits and
 * tracking statistics, as well as save/load, and not much else.
 *
 */

class BitStorage : public Storage
{
protected:
    std::vector<uint64_t> _tablesizes;
    size_t _n_tables;
    uint64_t _occupied_bins;
    uint64_t _n_unique_kmers;
    byte_t ** _counts;

public:
    BitStorage(const std::vector<uint64_t>& tablesizes) :
        _tablesizes(tablesizes)
    {
        _occupied_bins = 0;
        _n_unique_kmers = 0;

        _allocate_counters();
    }
    ~BitStorage()
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

    void _allocate_counters()
    {
        _n_tables = _tablesizes.size();

        _counts = new byte_t*[_n_tables];

        for (size_t i = 0; i < _n_tables; i++) {
            uint64_t tablesize = _tablesizes[i];
            uint64_t tablebytes = tablesize / 8 + 1;

            _counts[i] = new byte_t[tablebytes];
            memset(_counts[i], 0, tablebytes);
        }
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

    void save(std::string, uint16_t ksize);
    void load(std::string, uint16_t& ksize);

    // count number of occupied bins
    const uint64_t n_occupied() const
    {
        return _occupied_bins;
    }

    const uint64_t n_unique_kmers() const
    {
        return _n_unique_kmers;
    }

    // Get and set the hashbits for the given kmer hash.
    // Generally, it is better to keep tests and mutations separate,
    // but, in the interests of efficiency and thread safety,
    // tests and mutations are being blended here against conventional
    // software engineering wisdom.
    inline
    count_t 
    test_and_set_bits( hashing::hash_t khash )
    {
        bool is_new_kmer = false;

        for (size_t i = 0; i < _n_tables; i++) {
            uint64_t bin = khash % _tablesizes[i];
            uint64_t byte = bin / 8;
            unsigned char bit = (unsigned char)(1 << (bin % 8));

            unsigned char bits_orig = __sync_fetch_and_or( *(_counts + i) +
                                      byte, bit );
            if (!(bits_orig & bit)) {
                if (i == 0) {
                    __sync_add_and_fetch( &_occupied_bins, 1 );
                }
                is_new_kmer = true;
            }
        } // iteration over hashtables

        if (is_new_kmer) {
            __sync_add_and_fetch( &_n_unique_kmers, 1 );
            return 1; // kmer not seen before
        }

        return 0; // kmer already seen
    } // test_and_set_bits

    inline bool add(hashing::hash_t khash)
    {
        return test_and_set_bits(khash);
    }

    // get the count for the given k-mer hash.
    inline const count_t get_count(hashing::hash_t khash) const
    {
        for (size_t i = 0; i < _n_tables; i++) {
            uint64_t bin = khash % _tablesizes[i];
            uint64_t byte = bin / 8;
            unsigned char bit = bin % 8;

            if (!(_counts[i][byte] & (1 << bit))) {
                return 0;
            }
        }
        return 1;
    }

    // Writing to the tables outside of defined methods has undefined behavior!
    // As such, this should only be used to return read-only interfaces
    byte_t ** get_raw_tables()
    {
        return _counts;
    }

    void reset();

    void update_from(const BitStorage&);
};

}
}

#endif
