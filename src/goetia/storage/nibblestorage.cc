/**
 * (c) Camille Scott, 2019
 * File   : nibblestorage.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
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

#include "goetia/storage/nibblestorage.hh"

#include <cstring>
#include <errno.h>
#include <sstream> // IWYU pragma: keep
#include <fstream>
#include <iostream>

#include "goetia/goetia.hh"

using namespace std;
using namespace goetia;
using namespace goetia::storage;


const bool
NibbleStorage::insert(value_type khash)
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

const count_t
NibbleStorage::insert_and_query(value_type khash)
{
    if (insert(khash)) {
        return 1;
    }
    return query(khash);
}

// get the count for the given k-mer hash.
const count_t
NibbleStorage::query(value_type khash) const
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

void
NibbleStorage::save(std::string outfilename, uint16_t ksize)
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

    unsigned char ht_type = SAVED_SMALLCOUNT;
    outfile.write((const char *) &ht_type, 1);

    outfile.write((const char *) &save_ksize, sizeof(save_ksize));
    outfile.write((const char *) &save_n_tables, sizeof(save_n_tables));
    outfile.write((const char *) &save_occupied_bins,
                  sizeof(save_occupied_bins));

    for (unsigned int i = 0; i < save_n_tables; i++) {
        save_tablesize = _tablesizes[i];

        outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));
        outfile.write((const char *) _counts[i], save_tablesize / 2 + 1);
    }
}

void
NibbleStorage::load(std::string infilename, uint16_t& ksize)
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
            err = "Cannot open k-mer count file: " + infilename;
        } else {
            err = "Unknown error in opening file: " + infilename;
        }
        throw GoetiaFileException(err + " " + strerror(errno));
    } catch (const std::exception &e) {
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
        char signature [4];
        unsigned char version = 0, ht_type = 0;

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
                << " while reading k-mer count file from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw GoetiaFileException(err.str());
        } else if (!(ht_type == SAVED_SMALLCOUNT)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading k-mer count file from " << infilename;
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
            _counts[i] = NULL;
        }

        for (unsigned int i = 0; i < _n_tables; i++) {
            uint64_t tablesize;
            uint64_t tablebytes;

            infile.read((char *) &save_tablesize, sizeof(save_tablesize));

            tablebytes = save_tablesize / 2 + 1;
            tablesize = save_tablesize;
            _tablesizes.push_back(tablesize);

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
            err = "Unexpected end of k-mer count file: " + infilename;
        } else {
            err = "Error reading from k-mer count file: " + infilename + " "
                  + strerror(errno);
        }
        throw GoetiaFileException(err);
    } catch (const std::exception &e) {
        std::string err = "Error reading from k-mer count file: " + infilename + " "
                          + strerror(errno);
        throw GoetiaFileException(err);
    }
}

