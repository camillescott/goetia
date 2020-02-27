/**
 * (c) Camille Scott, 2019
 * File   : qfstorage.cc
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

#include "goetia/storage/qfstorage.hh"

#include <memory>
#include <errno.h>
#include <cstring>
#include <sstream> // IWYU pragma: keep
#include <fstream>
#include <iostream>

#include "goetia/goetia.hh"
#include "goetia/storage/cqf/gqf.h"

using namespace std;
using namespace goetia;
using namespace goetia::storage;


QFStorage::QFStorage(int size)
    : _size(size)
{
    cf = std::make_shared<QF>();
    // size is the power of two to specify the number of slots in
    // the filter (2**size). Third argument sets the number of bits used
    // in the key (current value of size+8 is copied from the CQF example)
    // Final argument is the number of bits allocated for the value, which
    // we do not use.
    qf_init(cf.get(), (1ULL << size), size+8, 0);
}


QFStorage::~QFStorage() 
{ 
    qf_destroy(cf.get());
}


std::shared_ptr<QFStorage>
QFStorage::clone() const {
    return std::make_shared<QFStorage>(_size);
}


const bool
QFStorage::insert(value_type khash) {
    bool is_new = query(khash) == 0;
    qf_insert(cf.get(), khash % cf->range, 0, 1);
    return is_new;
}


const count_t
QFStorage::insert_and_query(value_type khash) {
    qf_insert(cf.get(), khash % cf->range, 0, 1);
    return query(khash);
}


const count_t
QFStorage::query(value_type khash) const 
{
    return qf_count_key_value(cf.get(), khash % cf->range, 0);
}


std::vector<uint64_t>
QFStorage::get_tablesizes() const 
{ 
    return {cf->xnslots}; 
}


const uint64_t
QFStorage::n_unique_kmers() const 
{ 
    return cf->ndistinct_elts; 
}


const uint64_t
QFStorage::n_occupied() const 
{ 
    return cf->noccupied_slots; 
}


void
QFStorage::save(std::string outfilename, uint16_t ksize)
{
    ofstream outfile(outfilename.c_str(), ios::binary);

    unsigned char version = SAVED_FORMAT_VERSION;
    unsigned char ht_type = SAVED_QFCOUNT;

    outfile.write(SAVED_SIGNATURE, 4);
    outfile.write((const char *) &version, 1);
    outfile.write((const char *) &ht_type, 1);
    outfile.write((const char *) &ksize, sizeof(ksize));

    /* just a hack to handle __uint128_t value. Don't know a better to handle it
     * right now */
    uint64_t tmp_range;
    tmp_range = cf->range;

    outfile.write((const char *) &cf->nslots, sizeof(cf->nslots));
    outfile.write((const char *) &cf->xnslots, sizeof(cf->xnslots));
    outfile.write((const char *) &cf->key_bits, sizeof(cf->key_bits));
    outfile.write((const char *) &cf->value_bits, sizeof(cf->value_bits));
    outfile.write((const char *) &cf->key_remainder_bits, sizeof(cf->key_remainder_bits));
    outfile.write((const char *) &cf->bits_per_slot, sizeof(cf->bits_per_slot));
    outfile.write((const char *) &tmp_range, sizeof(tmp_range));
    outfile.write((const char *) &cf->nblocks, sizeof(cf->nblocks));
    outfile.write((const char *) &cf->nelts, sizeof(cf->nelts));
    outfile.write((const char *) &cf->ndistinct_elts, sizeof(cf->ndistinct_elts));
    outfile.write((const char *) &cf->noccupied_slots, sizeof(cf->noccupied_slots));

    #if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
        outfile.write((const char *) cf->blocks, sizeof(qfblock) * cf->nblocks);
    #else
        outfile.write((const char *) cf->blocks,
                      (sizeof(qfblock) + SLOTS_PER_BLOCK * cf->bits_per_slot / 8) * cf->nblocks);
    #endif
    outfile.close();
}


void
QFStorage::load(std::string infilename, uint16_t &ksize)
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
    uint16_t save_ksize = 0;
    char signature [4];
    unsigned char version = 0, ht_type = 0;
    uint64_t tmp_range;

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
    } else if (!(ht_type == SAVED_QFCOUNT)) {
        std::ostringstream err;
        err << "Incorrect file format type " << (int) ht_type
            << " expected " << (int) SAVED_QFCOUNT
            << " while reading k-mer count file from " << infilename;
        throw GoetiaFileException(err.str());
    }

    infile.read((char *) &save_ksize, sizeof(save_ksize));
    ksize = save_ksize;

    infile.read((char *) &cf->nslots, sizeof(cf->nslots));
    infile.read((char *) &cf->xnslots, sizeof(cf->xnslots));
    infile.read((char *) &cf->key_bits, sizeof(cf->key_bits));
    infile.read((char *) &cf->value_bits, sizeof(cf->value_bits));
    infile.read((char *) &cf->key_remainder_bits, sizeof(cf->key_remainder_bits));
    infile.read((char *) &cf->bits_per_slot, sizeof(cf->bits_per_slot));
    infile.read((char *) &tmp_range, sizeof(tmp_range));

    infile.read((char *) &cf->nblocks, sizeof(cf->nblocks));
    infile.read((char *) &cf->nelts, sizeof(cf->nelts));
    infile.read((char *) &cf->ndistinct_elts, sizeof(cf->ndistinct_elts));
    infile.read((char *) &cf->noccupied_slots, sizeof(cf->noccupied_slots));
    /* just a hack to handle __uint128_t value. Don't know a better to handle it
     * right now */
    cf->range = tmp_range;
    // deallocate previously allocated blocks
    free(cf->blocks);
    /* allocate the space for the actual qf blocks */
    #if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
        cf->blocks = (qfblock *)calloc(cf->nblocks, sizeof(qfblock));
    #else
        cf->blocks = (qfblock *)calloc(cf->nblocks, sizeof(qfblock) + SLOTS_PER_BLOCK * cf->bits_per_slot / 8);
    #endif
    #if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
        infile.read((char *) cf->blocks, sizeof(qfblock) * cf->nblocks);
    #else
        infile.read((char *) cf->blocks,
                    (sizeof(qfblock) + SLOTS_PER_BLOCK * cf->bits_per_slot / 8) * cf->nblocks);
    #endif
    infile.close();
}


