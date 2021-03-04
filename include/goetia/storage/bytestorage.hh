/**
 * (c) Camille Scott, 2019
 * File   : bytestorage.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 *
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

#ifndef GOETIA_BYTESTORAGE_HH
#define GOETIA_BYTESTORAGE_HH

#include <cassert>
#include <cstring>
#include <array>
#include <memory>
#include <mutex>
#include <unordered_map>

#include "goetia/meta.hh"
#include "goetia/storage/storage.hh"

#   define MAX_KCOUNT 255

namespace goetia {
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

class ByteStorage;
class ByteStorageFile;
class ByteStorageFileReader;
class ByteStorageFileWriter;
class ByteStorageGzFileReader;
class ByteStorageGzFileWriter;


template<>
struct StorageTraits<ByteStorage> {
    static constexpr bool is_probabilistic = true;
    static constexpr bool is_counting      = true;
 
    typedef std::tuple<uint64_t, uint16_t> params_type;
    static constexpr params_type default_params = std::make_tuple(1'000'000, 4);
};


class ByteStorage: public Storage<uint64_t>,
                   public Tagged<ByteStorage> {

    friend class ByteStorageFile;
    friend class ByteStorageFileReader;
    friend class ByteStorageFileWriter;
    friend class ByteStorageGzFileReader;
    friend class ByteStorageGzFileWriter;
    friend class CountGraph;
public:

    typedef uint64_t value_type;
    using Storage<uint64_t>::CountMap;
    using Traits = StorageTraits<ByteStorage>;

protected:

    count_t         _max_count;
    unsigned int    _max_bigcount;

    uint32_t _bigcount_spin_lock;
    std::vector<uint64_t> _tablesizes;
    size_t   _n_tables;
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
    CountMap _bigcounts;

    ByteStorage(uint64_t max_table, uint16_t N)
        : ByteStorage(get_n_primes_near_x(N, max_table))
    {
    }

    // constructor: create an empty CountMin sketch.
    ByteStorage(const std::vector<uint64_t>& tablesizes ) :
        _max_count(MAX_KCOUNT),
        _max_bigcount(MAX_BIGCOUNT),
        _bigcount_spin_lock(false), 
        _tablesizes(tablesizes),
        _n_unique_kmers(0), 
        _occupied_bins(0)
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

    std::shared_ptr<ByteStorage> clone() const {
        return std::make_shared<ByteStorage>(this->_tablesizes);
    }

    static std::shared_ptr<ByteStorage> build(uint64_t max_table, uint16_t N) {
        return std::make_shared<ByteStorage>(max_table, N);
    }

    static std::shared_ptr<ByteStorage> build(typename Traits::params_type params) {
        return make_shared_from_tuple<ByteStorage>(std::move(params));
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

    double estimated_fp() {
        double fp = n_occupied() / _tablesizes[0];
        fp = pow(fp, n_tables());
        return fp;
    }

    void save(std::string, uint16_t);
    void load(std::string, uint16_t&);

    const bool insert(value_type khash);

    const count_t insert_and_query(value_type khash);

    // get the count for the given k-mer hash.
    const count_t query(value_type khash) const;
    // Get direct access to the counts.
    //
    // Note:
    // Writing to the tables outside of defined methods has undefined behavior!
    // As such, this should only be used to return read-only interfaces
    byte_t ** get_raw_tables()
    {
        return _counts;
    }
    // not implemented
    static std::shared_ptr<ByteStorage> deserialize(std::ifstream& in) {
        return {};
    }

    void serialize(std::ofstream& out) {}
};


template<> 
struct is_probabilistic<ByteStorage> { 
      static const bool value = true;
};


template<>
struct is_counting<ByteStorage> {
    static const bool value = true;
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
