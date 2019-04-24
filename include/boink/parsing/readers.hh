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

#include <regex.h>
#include <stddef.h>
#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <memory>

#include "boink/boink.hh"
#include "boink/parsing/parsing.hh"


namespace seqan
{
    class SequenceStream; // forward dec seqan dep
}

namespace boink
{

namespace parsing
{

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

struct UnknownPairReadingMode : public  BoinkException {
    explicit UnknownPairReadingMode(const std::string& msg) :
        BoinkException(msg) {}
    UnknownPairReadingMode() :
        BoinkException("Unknown pair reading mode supplied.") {}
};

struct InvalidReadPair : public  BoinkException {
    explicit InvalidReadPair(const std::string& msg) :
        BoinkException(msg) {}
    InvalidReadPair() :
        BoinkException("Invalid read pair detected.") {}
};


template<typename SeqIO>
class ReadParser
{
protected:
    std::unique_ptr<SeqIO> _parser;
    regex_t _re_read_2_nosub;
    regex_t _re_read_1;
    regex_t _re_read_2;
    void _init();

    ReadPair _get_next_read_pair_in_ignore_mode();
    ReadPair _get_next_read_pair_in_error_mode();
    bool _is_valid_read_pair(
        ReadPair &the_read_pair,
        regmatch_t &match_1,
        regmatch_t &match_2
    );

public:
    enum {
        PAIR_MODE_IGNORE_UNPAIRED,
        PAIR_MODE_ERROR_ON_UNPAIRED
    };

    ReadParser(std::unique_ptr<SeqIO> pf);

    ReadParser(ReadParser& other);
    ReadParser& operator=(ReadParser& other);

    ReadParser(ReadParser&&) noexcept;
    ReadParser& operator=(ReadParser&&) noexcept;

    virtual ~ReadParser();

    Read get_next_read();
    ReadPair get_next_read_pair(uint8_t mode = PAIR_MODE_ERROR_ON_UNPAIRED);

    size_t get_num_reads();
    bool is_complete();
    void close();
}; // class ReadParser

// Alias for generic/templated ReadParser pointer
template<typename T> using ReadParserPtr = std::shared_ptr<ReadParser<T>>;
template<typename T> using WeakReadParserPtr = std::weak_ptr<ReadParser<T>>;

// Convenience function
template<typename SeqIO>
ReadParserPtr<SeqIO> get_parser(const std::string& filename);


class FastxReader
{
private:
    std::string _filename;
    std::unique_ptr<seqan::SequenceStream> _stream;
    uint32_t _spin_lock;
    size_t _num_reads;
    bool _have_qualities;
    void _init();

public:
    FastxReader();
    FastxReader(const std::string& infile);

    FastxReader(FastxReader& other);
    FastxReader& operator=(FastxReader& other);

    FastxReader(FastxReader&&) noexcept;
    FastxReader& operator=(FastxReader&&) noexcept;

    ~FastxReader();

    static std::shared_ptr<ReadParser<FastxReader>> build(const std::string& filename) {
        return get_parser<FastxReader>(filename);
    }

    Read get_next_read();
    bool is_complete();
    size_t get_num_reads();
    void close();
}; // class FastxReader


// Alias for instantiated ReadParsers
typedef std::shared_ptr<ReadParser<FastxReader>> FastxParserPtr;
typedef std::weak_ptr<ReadParser<FastxReader>> WeakFastxParserPtr;

/*
class BrokenPairedReader {

    std::unique_ptr<ReadParserPtr<FastxReader>> _parser;
    uint32_t                                    _min_length;
    bool                                        _force_single;
    bool                                        _require_paired;

    Read                                        _buffer_sequence;
    bool                                        _buffer_full;


    ReadBundle next() {
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


template <class ParserType = FastxReader>
class SplitPairedReader {

    ReadParserPtr<ParserType> left_parser;
    ReadParserPtr<ParserType> right_parser;
    uint32_t _min_length;
    bool     _force_name_match;
    uint64_t _n_reads;

public:

    SplitPairedReader(const std::string &left,
                 const std::string &right,
                 uint32_t min_length=0,
                 bool force_name_match=false)
        : _min_length(min_length),
          _force_name_match(force_name_match) {
        
        left_parser = get_parser<ParserType>(left);
        right_parser = get_parser<ParserType>(right);
    }

    bool is_complete() const {
        if (left_parser->is_complete() != right_parser->is_complete()) {
            throw BoinkException("Mismatched split paired files.");
        }
        return left_parser->is_complete();
    }

    ReadBundle next() {
        ReadBundle result;
        try {
            result.left = this->left_parser->get_next_read();
            result.right = this->right_parser->get_next_read();
        } catch (NoMoreReadsAvailable) {
            result.has_left = false;
            result.has_right = false;
            return result;
        }

        result.left.set_clean_seq();
        result.right.set_clean_seq();
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
