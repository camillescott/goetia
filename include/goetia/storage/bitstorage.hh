/**
 * (c) Camille Scott, 2019
 * File   : bitstorage.hh
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

#ifndef GOETIA_BITSTORAGE_HH
#define GOETIA_BITSTORAGE_HH

#include <cassert>
#include <cmath>
#include <cstring>
#include <array>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <utility>

#include "goetia/meta.hh"
#include "goetia/storage/storage.hh"


namespace goetia {
namespace storage {


class BitStorage;


template<>
struct StorageTraits<BitStorage> {
    static constexpr bool is_probabilistic = true;
    static constexpr bool is_counting      = false;

    typedef std::tuple<uint64_t, uint16_t> params_type;
    static constexpr params_type default_params = std::make_tuple(1'000'000, 4);
};


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

class BitStorage : public Storage<uint64_t>,
                   public Tagged<BitStorage>
{
protected:
    std::vector<uint64_t> _tablesizes;
    size_t   _n_tables;
    uint64_t _occupied_bins;
    uint64_t _n_unique_kmers;
    byte_t ** _counts;

public:

    using Storage<uint64_t>::value_type;
    using Traits = StorageTraits<BitStorage>;

    BitStorage(uint64_t max_table, uint16_t N)
        : BitStorage(get_n_primes_near_x(N, max_table))
    {
    }

    BitStorage(const std::vector<uint64_t>& tablesizes) :
        _tablesizes(tablesizes),
        _n_tables(tablesizes.size())
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

    std::shared_ptr<BitStorage> clone() const {
        return std::make_shared<BitStorage>(this->_tablesizes);
    }

    static std::shared_ptr<BitStorage> build(uint64_t max_table, uint16_t N);
    static std::shared_ptr<BitStorage> build(const typename Traits::params_type& params);

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
        return _tablesizes.size();
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

    double estimated_fp() {
        double fp = n_occupied() / _tablesizes[0];
        fp = pow(fp, n_tables());
        return fp;
    }

    // Get and set the hashbits for the given kmer hash.
    // Generally, it is better to keep tests and mutations separate,
    // but, in the interests of efficiency and thread safety,
    // tests and mutations are being blended here against conventional
    // software engineering wisdom.
    const bool insert( value_type khash );

    const count_t insert_and_query(value_type khash);

    // get the count for the given k-mer hash.
    const count_t query(value_type khash) const;

    // Writing to the tables outside of defined methods has undefined behavior!
    // As such, this should only be used to return read-only interfaces
    byte_t ** get_raw_tables()
    {
        return _counts;
    }

    void reset();

    void update_from(const BitStorage&);

    // not implemented
    static std::shared_ptr<BitStorage> deserialize(std::ifstream& in) {
        return {};
    }

    void serialize(std::ofstream& out) {}
};

template<>
struct is_probabilistic<BitStorage> { 
      static const bool value = true;
};

template<>
struct is_counting<BitStorage> {
    static const bool value = false;
};


}
}

#endif
