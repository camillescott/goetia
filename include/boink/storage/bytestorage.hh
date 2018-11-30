/* bytestorage.hh -- boink-modified oxli bytestorage
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

#ifndef BOINK_BYTESTORAGE_HH
#define BOINK_BYTESTORAGE_HH

#include <cassert>
#include <cstring>
#include <array>
#include <memory>
#include <mutex>
#include <unordered_map>

#include "boink/hashing/hashing_types.hh"
#include "boink/storage/storage.hh"

#   define MAX_KCOUNT 255

namespace boink {
namespace storage {

/*
 * \class ByteStorage
 *
 * \brief A CountMin sketch implementation.
 *
 * ByteStorage is used to track counts of k-mers by Hashtable
 * and derived classes.  It contains 'n_tables' different tables of
 * bytesizes specified in 'tablesizes'.
 *
 * Like other Storage classes, ByteStorage manages setting the bits and
 * tracking statistics, as well as save/load, and not much else.
 *
 */

class ByteStorageFile;
class ByteStorageFileReader;
class ByteStorageFileWriter;
class ByteStorageGzFileReader;
class ByteStorageGzFileWriter;

class ByteStorage : public Storage {
    friend class ByteStorageFile;
    friend class ByteStorageFileReader;
    friend class ByteStorageFileWriter;
    friend class ByteStorageGzFileReader;
    friend class ByteStorageGzFileWriter;
    friend class CountGraph;
protected:
    unsigned int    _max_count;
    unsigned int    _max_bigcount;

    uint32_t _bigcount_spin_lock;
    std::vector<uint64_t> _tablesizes;
    size_t _n_tables;
    uint64_t _n_unique_kmers;
    uint64_t _occupied_bins;

    byte_t ** _counts;

    // initialize counts with empty hashtables.
    void _allocate_counters()
    {
        _n_tables = _tablesizes.size();

        _counts = new byte_t*[_n_tables];
        for (size_t i = 0; i < _n_tables; i++) {
            _counts[i] = new byte_t[_tablesizes[i]];
            memset(_counts[i], 0, _tablesizes[i]);
        }
    }
public:
    KmerCountMap _bigcounts;

    // constructor: create an empty CountMin sketch.
    ByteStorage(const std::vector<uint64_t>& tablesizes ) :
        _max_count(MAX_KCOUNT), _max_bigcount(MAX_BIGCOUNT),
        _bigcount_spin_lock(false), _tablesizes(tablesizes),
        _n_unique_kmers(0), _occupied_bins(0)
    {
        _supports_bigcount = true;
        _allocate_counters();
    }

    // destructor: clear out the memory.
    ~ByteStorage()
    {
        if (_counts) {
            for (size_t i = 0; i < _n_tables; i++) {
                if (_counts[i]) {
                    delete[] _counts[i];
                    _counts[i] = NULL;
                }
            }

            delete[] _counts;
            _counts = NULL;

            _n_tables = 0;
        }
    }

    void reset()
    {
        for (unsigned int table_num = 0; table_num < _n_tables; table_num++) {
            uint64_t tablesize = _tablesizes[table_num];
            memset(_counts[table_num], 0, tablesize);
        }
    }

    std::vector<uint64_t> get_tablesizes() const
    {
        return _tablesizes;
    }

    const uint64_t n_unique_kmers() const
    {
        return _n_unique_kmers;
    }
    const size_t n_tables() const
    {
        return _n_tables;
    }
    const uint64_t n_occupied() const
    {
        return _occupied_bins;
    }

    void save(std::string, uint16_t);
    void load(std::string, uint16_t&);

    inline count_t test_and_set_bits(hashing::hash_t khash)
    {
        count_t x = get_count(khash);
        add(khash);
        return !x;
    }

    inline bool add(hashing::hash_t khash)
    {
        bool is_new_kmer = false;
        unsigned int  n_full	  = 0;

        // add one to each entry in each table.
        for (unsigned int i = 0; i < _n_tables; i++) {
            const uint64_t bin = khash % _tablesizes[i];
            byte_t current_count = _counts[ i ][ bin ];

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
            // NOTE: Technically, multiple threads can cause the bin to spill
            //	 over max_count a little, if they all read it as less than
            //	 max_count before any of them increment it.
            //	 However, do we actually care if there is a little
            //	 bit of slop here? It can always be trimmed off later, if
            //	 that would help with stats.

            if ( _max_count > current_count ) {
                __sync_add_and_fetch( *(_counts + i) + bin, 1 );
            } else {
                n_full++;
            }
        } // for each table

        // if all tables are full for this position, then add in bigcounts.
        if (n_full == _n_tables && _use_bigcount) {
            while (!__sync_bool_compare_and_swap(&_bigcount_spin_lock, 0, 1));
            if (_bigcounts[khash] == 0) {
                _bigcounts[khash] = _max_count + 1;
            } else {
                if (_bigcounts[khash] < _max_bigcount) {
                    _bigcounts[khash] += 1;
                }
            }
            __sync_bool_compare_and_swap( &_bigcount_spin_lock, 1, 0 );
        }

        if (is_new_kmer) {
            __sync_add_and_fetch(&_n_unique_kmers, 1);
        }

        return is_new_kmer;
    }

    // get the count for the given k-mer hash.
    inline const count_t get_count(hashing::hash_t khash) const
    {
        unsigned int	  max_count	= _max_count;
        count_t  min_count	= max_count; // bound count by max.

        // first, get the min count across all tables (standard CMS).
        for (unsigned int i = 0; i < _n_tables; i++) {
            count_t the_count = _counts[i][khash % _tablesizes[i]];
            if (the_count < min_count) {
                min_count = the_count;
            }
        }

        // if the count is saturated, check in the bigcount structure to
        // see if we've accumulated more counts.
        if (min_count == max_count && _use_bigcount) {
            KmerCountMap::const_iterator it = _bigcounts.find(khash);
            if (it != _bigcounts.end()) {
                min_count = it->second;
            }
        }
        return min_count;
    }
    // Get direct access to the counts.
    //
    // Note:
    // Writing to the tables outside of defined methods has undefined behavior!
    // As such, this should only be used to return read-only interfaces
    byte_t ** get_raw_tables()
    {
        return _counts;
    }

};

// Helper classes for saving ByteStorage objs to disk & loading them.

class ByteStorageFile
{
public:
    static void load(const std::string &infilename,
                     uint16_t &ksize,
                     ByteStorage &store);
    static void save(const std::string &outfilename,
                     const uint16_t ksize,
                     const ByteStorage &store);
};

class ByteStorageFileReader : public ByteStorageFile
{
public:
    ByteStorageFileReader(const std::string &infilename,
                          uint16_t &ksize,
                          ByteStorage &store);
};

class ByteStorageGzFileReader : public ByteStorageFile
{
public:
    ByteStorageGzFileReader(const std::string &infilename,
                            uint16_t &ksize,
                            ByteStorage &store);
};


class ByteStorageFileWriter : public ByteStorageFile
{
public:
    ByteStorageFileWriter(const std::string &outfilename,
                          const uint16_t ksize,
                          const ByteStorage &store);
};

class ByteStorageGzFileWriter : public ByteStorageFile
{
public:
    ByteStorageGzFileWriter(const std::string &outfilename,
                            const uint16_t ksize,
                            const ByteStorage &store);
};


}
}

#endif
