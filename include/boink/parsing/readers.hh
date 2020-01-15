/* readers.hh
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 * 
 **** END BOINK LICENSE BLOCK
 *
 * This file is part of khmer, https://github.com/dib-lab/khmer/, and is
 * Copyright (C) 2012-2015, Michigan State University.
 * Copyright (C) 2015-2016, The Regents of the University of California.
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
 *     * Neither the name of the Michigan State University nor the names
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

#ifndef BOINK_READERS_HH
#define BOINK_READERS_HH

#include <stddef.h>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <memory>
#include <zlib.h>

#include "boink/boink.hh"
#include "boink/parsing/parsing.hh"
#include "boink/sequences/alphabets.hh"

extern "C" {
#include "kseq.h"
}
KSEQ_INIT(gzFile, gzread)


namespace boink {
namespace parsing {

struct NoMoreReadsAvailable : public  BoinkFileException {
    explicit NoMoreReadsAvailable(const std::string& msg) :
        BoinkFileException(msg) {}
    NoMoreReadsAvailable() :
        BoinkFileException("No more reads available in this stream.") {}
};

struct InvalidRead : public  BoinkException {
    explicit InvalidRead(const std::string& msg) :
        BoinkException(msg) {}
    InvalidRead() :
        BoinkException("Invalid FASTA/Q read") {}
};

struct InvalidReadPair : public  BoinkException {
    explicit InvalidReadPair(const std::string& msg) :
        BoinkException(msg) {}
    InvalidReadPair() :
        BoinkException("Invalid read pair detected.") {}
};


template<typename ParserType>
class SequenceReader {

protected:

    std::unique_ptr<ParserType> _parser;

public:

    SequenceReader(std::unique_ptr<ParserType> pf);

    SequenceReader(SequenceReader& other);
    SequenceReader& operator=(SequenceReader& other);

    SequenceReader(SequenceReader&&) noexcept;
    SequenceReader& operator=(SequenceReader&&) noexcept;

    static std::shared_ptr<SequenceReader<ParserType>> build(const std::string& filename);

    virtual ~SequenceReader() {}

    Record next();

    size_t num_parsed() const;
    bool is_complete() const;
}; // class SequenceReader


template<class Alphabet = DNA_SIMPLE>
class FastxParser
{
private:
    std::string _filename;
    kseq_t *    _kseq;
    gzFile      _fp;
    uint32_t    _spin_lock;
    size_t      _num_parsed;
    bool        _have_qualities;
    bool        _is_complete;

public:
    FastxParser();
    FastxParser(const std::string& infile);

    FastxParser(FastxParser& other);
    FastxParser& operator=(FastxParser& other) = delete;;

    FastxParser(FastxParser&&) noexcept;
    FastxParser& operator=(FastxParser&&)  = delete;

    ~FastxParser();

    static std::shared_ptr<SequenceReader<FastxParser>> build(const std::string& filename) {
        return SequenceReader<FastxParser>::build(filename);
    }

    Record next();
    size_t num_parsed() const;
    bool is_complete() const;
}; // class FastxParser


/*
class BrokenPairedReader {

    std::unique_ptr<SequenceReaderPtr<FastxParser>> _parser;
    uint32_t                                    _min_length;
    bool                                        _force_single;
    bool                                        _require_paired;

    Read                                        _buffer_sequence;
    bool                                        _buffer_full;


    RecordBundle next() {
        Read first, second;
        bool is_pair;

        if (!_buffer_full) {
            try {
                first = this->parser->get_next_read();
            } catch (NoMoreReadsAvailable) {
                
            }

        } else {
            first = _buffer_sequence;
        }
    }

};
*/


template <class ParserType = FastxParser<>>
class SplitPairedReader {

    typedef SequenceReader<ParserType> reader_type;

    std::shared_ptr<reader_type> left_parser;
    std::shared_ptr<reader_type> right_parser;
    uint32_t                     _min_length;
    bool                         _force_name_match;
    uint64_t                     _n_reads;

public:

    SplitPairedReader(const std::string &left,
                      const std::string &right,
                      uint32_t min_length=0,
                      bool force_name_match=false)
        : _min_length(min_length),
          _force_name_match(force_name_match) {
        
        left_parser = reader_type::build(left);
        right_parser = reader_type::build(right);
    }

    bool is_complete() const {
        if (left_parser->is_complete() != right_parser->is_complete()) {
            throw BoinkException("Mismatched split paired files.");
        }
        return left_parser->is_complete();
    }

    RecordPair next() {
        RecordPair result;
        try {
            result.left = this->left_parser->next();
        } catch (NoMoreReadsAvailable) {
            result.has_left = false;
        }

        try {
            result.right = this->right_parser->next();
        } catch (NoMoreReadsAvailable) {
            result.has_right = false;
            return result;
        }

        result.has_left = true;
        result.has_right = true;

        if (this->_force_name_match) {
            if (!check_is_pair(result.left.name, result.right.name)) {
                throw BoinkException("Unpaired reads");
            }
        }

        if (this->_min_length > 0) {
            filter_length(result, this->_min_length);
        }

        return result;
    }
};


} // namespace parsing

} // namespace boink

#endif 
