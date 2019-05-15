/* bytestorage.cc -- boink-modified oxli storage
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

#include "boink/storage/bytestorage.hh"

#include <errno.h>
#include <sstream> // IWYU pragma: keep
#include <fstream>
#include <iostream>

#include "boink/boink.hh"
#include "boink/hashing/hashing_types.hh"
#include "zlib.h"

using namespace std;
using namespace boink;
using namespace boink::storage;
using namespace boink::hashing;


const bool
ByteStorage::insert(hashing::hash_t khash) {
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


const count_t
ByteStorage::query(hashing::hash_t khash) const
{
    count_t	 max_count	= _max_count;
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


const count_t
ByteStorage::insert_and_query(hashing::hash_t khash)
{
    if (insert(khash)) {
        // was new, return 1 from insert
        return 1;
    }
    return query(khash);
}


void ByteStorageFile::save(
    const std::string   &outfilename,
    uint16_t ksize,
    const ByteStorage &store)
{
    std::string filename(outfilename);
    size_t found = filename.find_last_of(".");
    std::string type = filename.substr(found + 1);

    if (type == "gz") {
        ByteStorageGzFileWriter(filename, ksize, store);
    } else {
        ByteStorageFileWriter(filename, ksize, store);
    }
}


ByteStorageFileReader::ByteStorageFileReader(
    const std::string   &infilename,
    uint16_t& ksize,
    ByteStorage &store)
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
        throw BoinkFileException(err + " " + strerror(errno));
    } catch (const std::exception &e) {
        std::string err = "Unknown error opening file: " + infilename + " "
                          + strerror(errno);
        throw BoinkFileException(err);
    }

    if (store._counts) {
        for (unsigned int i = 0; i < store._n_tables; i++) {
            delete[] store._counts[i];
            store._counts[i] = NULL;
        }
        delete[] store._counts;
        store._counts = NULL;
    }
    store._tablesizes.clear();

    try {
        unsigned int save_ksize = 0;
        unsigned char save_n_tables = 0;
        unsigned long long save_tablesize = 0;
        unsigned long long save_occupied_bins = 0;
        char signature [4];
        unsigned char version = 0, ht_type = 0, use_bigcount = 0;

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
            throw BoinkFileException(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading k-mer count file from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw BoinkFileException(err.str());
        } else if (!(ht_type == SAVED_COUNTING_HT)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading k-mer count file from " << infilename;
            throw BoinkFileException(err.str());
        }

        infile.read((char *) &use_bigcount, 1);
        infile.read((char *) &save_ksize, sizeof(save_ksize));
        infile.read((char *) &save_n_tables, sizeof(save_n_tables));
        infile.read((char *) &save_occupied_bins, sizeof(save_occupied_bins));

        ksize = (uint16_t) save_ksize;
        store._n_tables = (unsigned int) save_n_tables;
        store._occupied_bins = save_occupied_bins;

        store._use_bigcount = use_bigcount;

        store._counts = new byte_t*[store._n_tables];
        for (unsigned int i = 0; i < store._n_tables; i++) {
            store._counts[i] = NULL;
        }

        for (unsigned int i = 0; i < store._n_tables; i++) {
            uint64_t tablesize;

            infile.read((char *) &save_tablesize, sizeof(save_tablesize));

            tablesize = save_tablesize;
            store._tablesizes.push_back(tablesize);

            store._counts[i] = new byte_t[tablesize];

            unsigned long long loaded = 0;
            while (loaded != tablesize) {
                infile.read((char *) store._counts[i], tablesize - loaded);
                loaded += infile.gcount();
            }
        }

        uint64_t n_counts = 0;
        infile.read((char *) &n_counts, sizeof(n_counts));

        if (n_counts) {
            store._bigcounts.clear();

            hash_t kmer;
            count_t count;

            for (uint64_t n = 0; n < n_counts; n++) {
                infile.read((char *) &kmer, sizeof(kmer));
                infile.read((char *) &count, sizeof(count));
                store._bigcounts[kmer] = count;
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
        throw BoinkFileException(err);
    } catch (const std::exception &e) {
        std::string err = "Error reading from k-mer count file: " + infilename + " "
                          + strerror(errno);
        throw BoinkFileException(err);
    }
}

ByteStorageGzFileReader::ByteStorageGzFileReader(
    const std::string   &infilename,
    uint16_t &ksize,
    ByteStorage    &store)
{
    gzFile infile = gzopen(infilename.c_str(), "rb");
    if (infile == Z_NULL) {
        std::string err = "Cannot open k-mer count file: " + infilename;
        throw BoinkFileException(err);
    }

    if (store._counts) {
        for (unsigned int i = 0; i < store._n_tables; i++) {
            delete[] store._counts[i];
            store._counts[i] = NULL;
        }
        delete[] store._counts;
        store._counts = NULL;
    }
    store._tablesizes.clear();

    unsigned int save_ksize = 0;
    unsigned char save_n_tables = 0;
    unsigned long long save_tablesize = 0;
    unsigned long long save_occupied_bins = 0;
    char signature [4];
    unsigned char version, ht_type, use_bigcount;

    int read_s = gzread(infile, signature, 4);
    int read_v = gzread(infile, (char *) &version, 1);
    int read_t = gzread(infile, (char *) &ht_type, 1);
    if (read_s <= 0 || read_v <= 0 || read_t <= 0) {
        std::string err = "K-mer count file read error: " + infilename + " "
                          + strerror(errno);
        gzclose(infile);
        throw BoinkFileException(err);
    } else if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
        std::ostringstream err;
        err << "Does not start with signature for a oxli " <<
            "file: " << signature << " Should be: " <<
            SAVED_SIGNATURE;
        throw BoinkFileException(err.str());
    } else if (!(version == SAVED_FORMAT_VERSION)
               || !(ht_type == SAVED_COUNTING_HT)) {
        if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading k-mer count file from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            gzclose(infile);
            throw BoinkFileException(err.str());
        } else if (!(ht_type == SAVED_COUNTING_HT)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading k-mer count file from " << infilename;
            gzclose(infile);
            throw BoinkFileException(err.str());
        }
    }

    int read_b = gzread(infile, (char *) &use_bigcount, 1);
    int read_k = gzread(infile, (char *) &save_ksize, sizeof(save_ksize));
    int read_nt = gzread(infile, (char *) &save_n_tables,
                         sizeof(save_n_tables));
    int read_ob = gzread(infile, (char *) &save_occupied_bins,
                         sizeof(save_occupied_bins));

    if (read_b <= 0 || read_k <= 0 || read_nt <= 0 || read_ob <= 0) {
        std::string err = "K-mer count file header read error: " + infilename
                          + " " + strerror(errno);
        gzclose(infile);
        throw BoinkFileException(err);
    }

    ksize = (uint16_t) save_ksize;
    store._occupied_bins = save_occupied_bins;
    store._n_tables = (unsigned int) save_n_tables;

    store._use_bigcount = use_bigcount;

    store._counts = new byte_t*[store._n_tables];
    for (unsigned int i = 0; i < store._n_tables; i++) {
        uint64_t tablesize;

        read_b = gzread(infile, (char *) &save_tablesize,
                        sizeof(save_tablesize));

        if (read_b <= 0) {
            std::string gzerr = gzerror(infile, &read_b);
            std::string err = "K-mer count file header read error: "
                              + infilename;
            if (read_b == Z_ERRNO) {
                err = err + " " + strerror(errno);
            } else {
                err = err + " " + gzerr;
            }
            gzclose(infile);
            throw BoinkFileException(err);
        }

        tablesize = save_tablesize;
        store._tablesizes.push_back(tablesize);

        store._counts[i] = new byte_t[tablesize];

        uint64_t loaded = 0;
        while (loaded != tablesize) {
            unsigned long long  to_read_ll = tablesize - loaded;
            unsigned int        to_read_int;
            // Zlib can only read chunks of at most INT_MAX bytes.
            if (to_read_ll > INT_MAX) {
                to_read_int = INT_MAX;
            } else {
                to_read_int = (unsigned int) to_read_ll;
            }
            read_b = gzread(infile, (char *) store._counts[i], to_read_int);

            if (read_b <= 0) {
                std::string gzerr = gzerror(infile, &read_b);
                std::string err = "K-mer count file read error: " + infilename;
                if (read_b == Z_ERRNO) {
                    err = err + " " + strerror(errno);
                } else {
                    err = err + " " + gzerr;
                }
                gzclose(infile);
                throw BoinkFileException(err);
            }

            loaded += read_b;
        }
    }

    uint64_t n_counts = 0;
    read_b = gzread(infile, (char *) &n_counts, sizeof(n_counts));
    if (read_b <= 0) {
        std::string gzerr = gzerror(infile, &read_b);
        std::string err = "K-mer count header read error: " + infilename;
        if (read_b == Z_ERRNO) {
            err = err + " " + strerror(errno);
        } else {
            err = err + " " + gzerr;
        }
        gzclose(infile);
        throw BoinkFileException(err);
    }

    if (n_counts) {
        store._bigcounts.clear();

        hash_t kmer;
        count_t count;

        for (uint64_t n = 0; n < n_counts; n++) {
            int read_k = gzread(infile, (char *) &kmer, sizeof(kmer));
            int read_c = gzread(infile, (char *) &count, sizeof(count));

            if (read_k <= 0 || read_c <= 0) {
                std::string gzerr = gzerror(infile, &read_b);
                std::string err = "K-mer count read error: " + infilename;
                if (read_b == Z_ERRNO) {
                    err = err + " " + strerror(errno);
                } else {
                    err = err + " " + gzerr;
                }
                gzclose(infile);
                throw BoinkFileException(err);
            }

            store._bigcounts[kmer] = count;
        }
    }

    gzclose(infile);
}

ByteStorageFileWriter::ByteStorageFileWriter(
    const std::string   &outfilename,
    const uint16_t ksize,
    const ByteStorage& store)
{
    if (!store._counts[0]) {
        throw BoinkException();
    }

    unsigned int save_ksize = ksize;
    unsigned char save_n_tables = store._n_tables;
    unsigned long long save_tablesize;
    unsigned long long save_occupied_bins = store._occupied_bins;

    ofstream outfile(outfilename.c_str(), ios::binary);

    outfile.write(SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_COUNTING_HT;
    outfile.write((const char *) &ht_type, 1);

    unsigned char use_bigcount = 0;
    if (store._use_bigcount) {
        use_bigcount = 1;
    }
    outfile.write((const char *) &use_bigcount, 1);

    outfile.write((const char *) &save_ksize, sizeof(save_ksize));
    outfile.write((const char *) &save_n_tables, sizeof(save_n_tables));
    outfile.write((const char *) &save_occupied_bins,
                  sizeof(save_occupied_bins));

    for (unsigned int i = 0; i < save_n_tables; i++) {
        save_tablesize = store._tablesizes[i];

        outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));
        outfile.write((const char *) store._counts[i], save_tablesize);
    }

    uint64_t n_counts = store._bigcounts.size();
    outfile.write((const char *) &n_counts, sizeof(n_counts));

    if (n_counts) {
        KmerCountMap::const_iterator it = store._bigcounts.begin();

        for (; it != store._bigcounts.end(); ++it) {
            outfile.write((const char *) &it->first, sizeof(it->first));
            outfile.write((const char *) &it->second, sizeof(it->second));
        }
    }
    if (outfile.fail()) {
        throw BoinkFileException(strerror(errno));
    }
    outfile.close();
}

ByteStorageGzFileWriter::ByteStorageGzFileWriter(
    const std::string   &outfilename,
    const uint16_t ksize,
    const ByteStorage &store)
{
    if (!store._counts[0]) {
        throw BoinkException();
    }

    int errnum = 0;
    unsigned int save_ksize = ksize;
    unsigned char save_n_tables = store._n_tables;
    unsigned long long save_tablesize;
    unsigned long long save_occupied_bins = store._occupied_bins;

    gzFile outfile = gzopen(outfilename.c_str(), "wb");
    if (outfile == NULL) {
        const char * error = gzerror(outfile, &errnum);
        if (errnum == Z_ERRNO) {
            throw BoinkFileException(strerror(errno));
        } else {
            throw BoinkFileException(error);
        }
    }

    gzwrite(outfile, SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    gzwrite(outfile, (const char *) &version, 1);

    unsigned char ht_type = SAVED_COUNTING_HT;
    gzwrite(outfile, (const char *) &ht_type, 1);

    unsigned char use_bigcount = 0;
    if (store._use_bigcount) {
        use_bigcount = 1;
    }
    gzwrite(outfile, (const char *) &use_bigcount, 1);

    gzwrite(outfile, (const char *) &save_ksize, sizeof(save_ksize));
    gzwrite(outfile, (const char *) &save_n_tables, sizeof(save_n_tables));
    gzwrite(outfile, (const char *) &save_occupied_bins,
            sizeof(save_occupied_bins));

    for (unsigned int i = 0; i < save_n_tables; i++) {
        save_tablesize = store._tablesizes[i];

        gzwrite(outfile, (const char *) &save_tablesize,
                sizeof(save_tablesize));
        unsigned long long written = 0;
        while (written != save_tablesize) {
            unsigned long long  to_write_ll = save_tablesize - written;
            unsigned int        to_write_int;
            int                 gz_result;
            // Zlib can only write chunks of at most INT_MAX bytes.
            if (to_write_ll > INT_MAX) {
                to_write_int = INT_MAX;
            } else {
                to_write_int = (unsigned int) to_write_ll;
            }
            gz_result = gzwrite(outfile, (const char *) store._counts[i],
                                to_write_int);
            // Zlib returns 0 on error
            if (gz_result == 0) {
                int errcode = 0;
                const char *err_msg;
                std::ostringstream msg;

                msg << "gzwrite failed while writing counting hash: ";
                // Get zlib error
                err_msg = gzerror(outfile, &errcode);
                if (errcode != Z_ERRNO) {
                    // Zlib error, not stdlib
                    msg << err_msg;
                    gzclearerr(outfile);
                } else {
                    // stdlib error
                    msg << strerror(errno);
                }
                gzclose(outfile);
                throw BoinkFileException(msg.str());
            }
            written += gz_result;
        }
    }

    uint64_t n_counts = store._bigcounts.size();
    gzwrite(outfile, (const char *) &n_counts, sizeof(n_counts));

    if (n_counts) {
        KmerCountMap::const_iterator it = store._bigcounts.begin();

        for (; it != store._bigcounts.end(); ++it) {
            gzwrite(outfile, (const char *) &it->first, sizeof(it->first));
            gzwrite(outfile, (const char *) &it->second, sizeof(it->second));
        }
    }
    const char * error = gzerror(outfile, &errnum);
    if (errnum == Z_ERRNO) {
        throw BoinkFileException(strerror(errno));
    } else if (errnum != Z_OK) {
        throw BoinkFileException(error);
    }
    gzclose(outfile);
}

void ByteStorageFile::load(
    const std::string   &infilename,
    uint16_t &ksize,
    ByteStorage &store)
{
    std::string filename(infilename);
    size_t found = filename.find_last_of(".");
    std::string type = filename.substr(found + 1);

    if (type == "gz") {
        ByteStorageGzFileReader(filename, ksize, store);
    } else {
        ByteStorageFileReader(filename, ksize, store);
    }
}

void ByteStorage::save(std::string outfilename, uint16_t ksize)
{
    ByteStorageFile::save(outfilename, ksize, *this);
}

void ByteStorage::load(std::string infilename, uint16_t& ksize)
{
    ByteStorageFile::load(infilename, ksize, *this);
}

