/**
 * (c) Camille Scott, 2019
 * File   : bitstorage.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */
/* bitstorage.cc -- goetia-modified oxli storage
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
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

#include "goetia/storage/bitstorage.hh"

#include <errno.h>
#include <sstream> // IWYU pragma: keep
#include <fstream>
#include <iostream>

#include "goetia/goetia.hh"
#include "zlib.h"

using namespace std;
using namespace goetia;
using namespace goetia::storage;

std::shared_ptr<BitStorage>
BitStorage::build(uint64_t max_table, uint16_t N) {
    return std::make_shared<BitStorage>(max_table, N);
}


std::shared_ptr<BitStorage>
BitStorage::build(const typename StorageTraits<BitStorage>::params_type& params) {
    return make_shared_from_tuple<BitStorage>(params);
}


const bool
BitStorage::insert( value_type khash ) {
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


const count_t
BitStorage::insert_and_query(value_type khash)
{
    insert(khash);
    // presence filter, should always be 1 after insert
    return 1;
}

// get the count for the given k-mer hash.
const count_t
BitStorage::query(value_type khash) const
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


void
BitStorage::update_from(const BitStorage& other)
{
    if (_tablesizes != other._tablesizes) {
        throw GoetiaException("both nodegraphs must have same table sizes");
    }

    byte_t tmp = 0;
    for (unsigned int table_num = 0; table_num < _n_tables; table_num++) {
        byte_t * me = _counts[table_num];
        byte_t * ot = other._counts[table_num];
        uint64_t tablesize = _tablesizes[table_num];
        uint64_t tablebytes = tablesize / 8 + 1;

        for (uint64_t index = 0; index < tablebytes; index++) {
            // Bloom filters can be unioned with bitwise OR.
            // First, get the new value
            tmp = me[index] | ot[index];
            if (table_num == 0) {
                // We'd like for the merged filter to have an accurate
                // count of occupied bins.  First, observe that
                // HammingDistance(x,y) is equivalent to
                // HammingWeight(x^y).  Then, observe that the number
                // of additional occupied bins from the update is the
                // hamming distance between the original bin and the
                // OR'd bin. Thus, we can use the builtin popcountll
                // function, which calls a hardware instruction for
                // hamming weight, with the original and merged bin,
                // to find the number of additional occupied bins.
                _occupied_bins += __builtin_popcountll(me[index] ^ tmp);
            }
            me[index] = tmp;
        }
    }
}


void
BitStorage::reset()
{
    for (unsigned int table_num = 0; table_num < _n_tables; table_num++) {
        uint64_t tablesize = _tablesizes[table_num];
        uint64_t tablebytes = tablesize / 8 + 1;
        memset(_counts[table_num], 0, tablebytes);
    }
}


void
BitStorage::save(std::string outfilename, uint16_t ksize)
{
    if (!_counts[0]) {
        throw GoetiaException();
    }

    unsigned int save_ksize = ksize;
    unsigned char save_n_tables = _n_tables;
    unsigned long long save_tablesize;
    unsigned long long save_occupied_bins = _occupied_bins;

    ofstream outfile(outfilename.c_str(), ios::binary);

    outfile.write(SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_HASHBITS;
    outfile.write((const char *) &ht_type, 1);

    outfile.write((const char *) &save_ksize, sizeof(save_ksize));
    outfile.write((const char *) &save_n_tables, sizeof(save_n_tables));
    outfile.write((const char *) &save_occupied_bins,
                  sizeof(save_occupied_bins));

    for (unsigned int i = 0; i < _n_tables; i++) {
        save_tablesize = _tablesizes[i];
        unsigned long long tablebytes = save_tablesize / 8 + 1;

        outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));

        outfile.write((const char *) _counts[i], tablebytes);
    }
    if (outfile.fail()) {
        throw GoetiaFileException(strerror(errno));
    }
    outfile.close();
}

/**
 * Loads @param infilename into BitStorage, with error checking on
 * file type and file version.  Populates _counts internally.
 */
void
BitStorage::load(std::string infilename, uint16_t &ksize)
{
    ifstream infile;

    // configure ifstream to raise exceptions for everything.
    infile.exceptions(std::ifstream::failbit | std::ifstream::badbit |
                      std::ifstream::eofbit);

    try {
        infile.open(infilename.c_str(), ios::binary);
    } catch (std::ifstream::failure &e) {
        std::string err;
        if (!infile.is_open()) {
            err = "Cannot open k-mer graph file: " + infilename;
        } else {
            err = "Unknown error in opening file: " + infilename;
        }
        throw GoetiaFileException(err);
    } catch (const std::exception &e) {
        // Catching std::exception is a stopgap for
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66145
        std::string err = "Unknown error opening file: " + infilename + " "
                          + strerror(errno);
        throw GoetiaFileException(err);
    }

    if (_counts) {
        for (unsigned int i = 0; i < _n_tables; i++) {
            delete[] _counts[i];
            _counts[i] = NULL;
        }
        delete[] _counts;
        _counts = NULL;
    }
    _tablesizes.clear();

    try {
        unsigned int save_ksize = 0;
        unsigned char save_n_tables = 0;
        unsigned long long save_tablesize = 0;
        unsigned long long save_occupied_bins = 0;
        char signature[4];
        unsigned char version, ht_type;

        infile.read(signature, 4);
        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
            std::ostringstream err;
            err << "Does not start with signature for a oxli file: 0x";
            for(size_t i=0; i < 4; ++i) {
                err << std::hex << (int) signature[i];
            }
            err << " Should be: " << SAVED_SIGNATURE;
            throw GoetiaFileException(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading k-mer graph from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw GoetiaFileException(err.str());
        } else if (!(ht_type == SAVED_HASHBITS)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading k-mer graph from " << infilename;
            throw GoetiaFileException(err.str());
        }

        infile.read((char *) &save_ksize, sizeof(save_ksize));
        infile.read((char *) &save_n_tables, sizeof(save_n_tables));
        infile.read((char *) &save_occupied_bins, sizeof(save_occupied_bins));

        ksize = (uint16_t) save_ksize;
        _n_tables = (unsigned int) save_n_tables;
        _occupied_bins = save_occupied_bins;

        _counts = new byte_t*[_n_tables];
        for (unsigned int i = 0; i < _n_tables; i++) {
            uint64_t tablesize;
            unsigned long long tablebytes;

            infile.read((char *) &save_tablesize, sizeof(save_tablesize));

            tablesize = save_tablesize;
            _tablesizes.push_back(tablesize);

            tablebytes = tablesize / 8 + 1;
            _counts[i] = new byte_t[tablebytes];

            unsigned long long loaded = 0;
            while (loaded != tablebytes) {
                infile.read((char *) _counts[i], tablebytes - loaded);
                loaded += infile.gcount();
            }
        }
        infile.close();
    } catch (std::ifstream::failure &e) {
        std::string err;
        if (infile.eof()) {
            err = "Unexpected end of k-mer graph file: " + infilename;
        } else {
            err = "Error reading from k-mer graph file: " + infilename;
        }
        throw GoetiaFileException(err);
    } catch (const std::exception &e) {
        // Catching std::exception is a stopgap for
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66145
        std::string err = "Unknown error opening file: " + infilename + " "
                          + strerror(errno);
        throw GoetiaFileException(err);
    }
}
