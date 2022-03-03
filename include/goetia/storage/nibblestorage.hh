/**
 * (c) Camille Scott, 2019
 * File   : nibblestorage.hh
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

#ifndef GOETIA_NIBBLESTORAGE_HH
#define GOETIA_NIBBLESTORAGE_HH

#include <cassert>
#include <cstring>
#include <array>
#include <memory>
#include <mutex>
#include <unordered_map>

#include "goetia/meta.hh"
#include "goetia/storage/storage.hh"


namespace goetia {


class NibbleStorage;

template<>
struct StorageTraits<NibbleStorage> {
    static constexpr bool is_probabilistic = true;
    static constexpr bool is_counting      = true;
    static constexpr int  bits_per_slot    = 4;

    typedef std::tuple<uint64_t, uint16_t> params_type;
    static constexpr params_type default_params = std::make_tuple(1'000'000, 4);
};


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
class NibbleStorage : public Storage<uint64_t>,
                      public Tagged<NibbleStorage>
{
public:

    using Storage<uint64_t>::value_type;
    using Traits = StorageTraits<NibbleStorage>;

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
    uint64_t _table_index(const value_type k, const uint64_t tablesize) const
    {
        return (k % tablesize) / 2;
    }
    // Compute which half of the byte to use for this hash value
    uint8_t _mask(const value_type k, const uint64_t tablesize) const
    {
        return (k%tablesize)%2 ? 15 : 240;
    }
    // Compute which half of the byte to use for this hash value
    uint8_t _shift(const value_type k, const uint64_t tablesize) const
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

    static std::shared_ptr<NibbleStorage> build(uint64_t max_table, uint16_t N) {
        return std::make_shared<NibbleStorage>(max_table, N);
    }

    static std::shared_ptr<NibbleStorage> build(typename Traits::params_type params) {
        return make_shared_from_tuple<NibbleStorage>(std::move(params));
    }

    std::shared_ptr<NibbleStorage> clone() const {
        return std::make_shared<NibbleStorage>(this->_tablesizes);
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

    const bool insert(value_type khash);

    const count_t insert_and_query(value_type khash);

    // get the count for the given k-mer hash.
    const count_t query(value_type khash) const;

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
        double fp = static_cast<double>(n_occupied()) / 
                    static_cast<double>(_tablesizes[0]);
        fp = pow(fp, n_tables());
        return fp;
    }
    void save(std::string outfilename, uint16_t ksize);
    void load(std::string infilename, uint16_t& ksize);

    byte_t ** get_raw_tables()
    {
        return _counts;
    }

    // not implemented
    static std::shared_ptr<NibbleStorage> deserialize(std::ifstream& in) {
        return {};
    }

    void serialize(std::ofstream& out) {}
};

}

#endif
