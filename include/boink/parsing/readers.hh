/**
 * (c) Camille Scott, 2019
 * File   : readers.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.01.2020
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
#include <fstream>
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

    typedef Record value_type;

    SequenceReader(std::unique_ptr<ParserType> pf)
    {
        _parser = std::move(pf);
    }

    SequenceReader(SequenceReader& other)
    {
        _parser = std::move(other._parser);
    }

    SequenceReader& operator=(SequenceReader& other) {
        _parser = std::move(other._parser);
        return *this;
    }

    SequenceReader(SequenceReader&&) noexcept {}
    SequenceReader& operator=(SequenceReader&&) noexcept {}

    static auto build(const std::string& filename)
    -> std::shared_ptr<SequenceReader<ParserType>> {
        return std::make_shared<SequenceReader<ParserType>>(
                   std::move(std::make_unique<ParserType>(filename))
               );
    }

    virtual ~SequenceReader() {}

    Record next() {
        return _parser->next();
    }

    size_t num_parsed() const {
        return _parser->num_parsed();
    }

    bool is_complete() const {
        return _parser->is_complete();
    }
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

    typedef Record   value_type;
    typedef Alphabet alphabet;

    FastxParser() 
        : FastxParser("-")
    {
    }

    FastxParser(const std::string& infile)
        : _filename(infile),
          _spin_lock(0),
          _num_parsed(0),
          _have_qualities(false),
          _is_complete(false)
    {
        _fp = gzopen(_filename.c_str(), "r");
        _kseq = kseq_init(_fp);

        __asm__ __volatile__ ("" ::: "memory");
    }

    FastxParser(FastxParser& other)
        : _filename(std::move(other._filename)),
          _spin_lock(other._spin_lock),
          _num_parsed(other._num_parsed),
          _have_qualities(other._have_qualities),
          _fp(other._fp),
          _kseq(other._kseq),
          _is_complete(other._is_complete)
    {
        other._is_complete = true;
    }

    FastxParser& operator=(FastxParser& other) = delete;

    FastxParser(FastxParser&&) noexcept;
    FastxParser& operator=(FastxParser&&)  = delete;

    ~FastxParser() {
        kseq_destroy(_kseq);
        gzclose(_fp);
    }

    static std::shared_ptr<SequenceReader<FastxParser>> build(const std::string& filename) {
        return SequenceReader<FastxParser>::build(filename);
    }

    Record next() {
        Record record;
        
        //while (!__sync_bool_compare_and_swap(&_spin_lock, 0, 1));

        int stat = kseq_read(_kseq);

        if (stat >= 0) {
            try {
                Alphabet::validate(_kseq->seq.s, _kseq->seq.l);
            } catch (InvalidCharacterException &e) {
                _spin_lock = 0;
                throw e;
            } catch (std::exception& e) {
                _spin_lock = 0;
                throw e;
            }

            record.sequence.assign(_kseq->seq.s, _kseq->seq.l);
            record.name.assign(_kseq->name.s, _kseq->name.l);
            if (_kseq->qual.l) {
                record.quality.assign(_kseq->qual.s, _kseq->qual.l);
                if (_num_parsed == 0) {
                    _have_qualities = true;
                }
            }
            _num_parsed++;
        }

        //__asm__ __volatile__ ("" ::: "memory");
        //_spin_lock = 0;

        if (stat == -1) {
            _is_complete = true;
            throw NoMoreReadsAvailable();
        }

        if (stat == -2) {
            throw InvalidRead("Sequence and quality lengths differ");
        }

        if (stat == -3) {
            throw BoinkFileException("Error reading stream.");
        }

        return record;
    }

    size_t num_parsed() const {
        return _num_parsed;
    }

    bool is_complete() const {
        return _is_complete;
    }
}; // class FastxParser


template <class ParserType = FastxParser<>>
class SplitPairedReader {

    typedef SequenceReader<ParserType> reader_type;

    std::shared_ptr<reader_type> left_parser;
    std::shared_ptr<reader_type> right_parser;
    uint32_t                     _min_length;
    bool                         _force_name_match;
    uint64_t                     _n_reads;

public:

    typedef RecordPair value_type;

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

        std::exception_ptr left_exc_ptr;
        std::exception_ptr right_exc_ptr;

        try {
            result.left = this->left_parser->next();
        } catch (NoMoreReadsAvailable) {
            result.has_left = false;
        } catch (InvalidCharacterException) {
            result.has_left = false;
            left_exc_ptr = std::current_exception();
        } catch (InvalidRead) {
            result.has_left = false;
            left_exc_ptr = std::current_exception();
        }

        try {
            result.right = this->right_parser->next();
        } catch (NoMoreReadsAvailable) {
            result.has_right = false;
        } catch (InvalidCharacterException) {
            result.has_right = false;
            right_exc_ptr = std::current_exception();
        } catch (InvalidRead) {
            result.has_right = false;
            right_exc_ptr = std::current_exception();
        }

        if (!(result.has_left && result.has_right)) {
            if (left_exc_ptr) {
                std::rethrow_exception(left_exc_ptr);
            }
            if (right_exc_ptr) {
                std::rethrow_exception(right_exc_ptr);
            }
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

extern template class parsing::FastxParser<DNA_SIMPLE>;
extern template class parsing::FastxParser<DNAN_SIMPLE>;
extern template class parsing::SequenceReader<parsing::FastxParser<DNA_SIMPLE>>;
extern template class parsing::SequenceReader<parsing::FastxParser<DNAN_SIMPLE>>;
extern template class parsing::SplitPairedReader<parsing::FastxParser<DNA_SIMPLE>>;
extern template class parsing::SplitPairedReader<parsing::FastxParser<DNAN_SIMPLE>>;


} // namespace parsing

} // namespace boink

#endif 
