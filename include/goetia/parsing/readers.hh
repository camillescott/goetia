/**
 * (c) Camille Scott, 2019
 * File   : readers.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.01.2020
 *
 **** END GOETIA LICENSE BLOCK
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

#ifndef GOETIA_READERS_HH
#define GOETIA_READERS_HH

#include <stddef.h>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <memory>
#include <zlib.h>

#include "goetia/goetia.hh"
#include "goetia/parsing/parsing.hh"
#include "goetia/sequences/alphabets.hh"
#include "goetia/errors.hh"


extern "C" {
#include "kseq.h"
}
KSEQ_INIT(gzFile, gzread)


namespace goetia {


template<class Alphabet = DNA_SIMPLE>
class FastxParser
{
private:
    std::string _filename;
    kseq_t *    _kseq;
    gzFile      _fp;
    uint32_t    _spin_lock;
    size_t      _n_parsed;
    bool        _have_qualities;
    bool        _is_complete;
    bool        _strict;
    uint64_t    _n_skipped;

    uint32_t    _min_length;

public:

    typedef Record   value_type;
    typedef Alphabet alphabet;
    
    FastxParser() 
        : FastxParser("-")
    {
    }

    FastxParser(const std::string& infile,
               bool strict = false,
               uint32_t min_length = 0);


    FastxParser(FastxParser&& other)
        : _filename(std::move(other._filename)),
          _spin_lock(other._spin_lock),
          _n_parsed(other._n_parsed),
          _have_qualities(other._have_qualities),
          _fp(other._fp),
          _kseq(other._kseq),
          _is_complete(other._is_complete),
          _strict(other._strict),
          _n_skipped(other._n_skipped),
          _min_length(other._min_length)
    {
        other._is_complete = true;
    }

    FastxParser& operator=(FastxParser& other) = delete;
    ~FastxParser();

    static std::shared_ptr<FastxParser> build(const std::string& filename,
                                              bool strict = false,
                                              uint32_t min_length = 0) {
        return std::make_shared<FastxParser>(filename, strict, min_length);
    }

    std::optional<Record> next() {
        if (is_complete()) {
            throw std::out_of_range("No more reads available.");
        }

        Record record;
        
        //while (!__sync_bool_compare_and_swap(&_spin_lock, 0, 1));

        int stat = kseq_read(_kseq);

        if (stat >= 0) {

            if (!Alphabet::sanitize(_kseq->seq.s, _kseq->seq.l)) {
                _spin_lock = 0;
                ++_n_skipped;
                if (_strict) {
                    throw InvalidSequence(_n_parsed, _kseq->seq.s);
                } else {
                    stat = -4;
                }
            }

            if (_kseq->seq.l < _min_length) {
                stat = -4;
                ++_n_skipped;
            }

            if (stat >= 0) {
                record.sequence.assign(_kseq->seq.s, _kseq->seq.l);
                record.name.assign(_kseq->name.s, _kseq->name.l);
                if (_kseq->qual.l) {
                    record.quality.assign(_kseq->qual.s, _kseq->qual.l);
                    if (_n_parsed == 0) {
                        _have_qualities = true;
                    }
                }
            }
            ++_n_parsed;
        }

        //__asm__ __volatile__ ("" ::: "memory");
        //_spin_lock = 0;

        if (stat == -1) {
            _is_complete = true;
            return {};
        }

        if (stat == -2) {
            ++_n_skipped;
            throw InvalidRecord(_n_parsed, _kseq->seq.s);
        }

        if (stat == -3) {
            throw StreamReadError(_n_parsed, _filename);
        }

        if (stat == -4) {
            return {};
        }

        return record;
    }

    size_t n_parsed() const {
        return _n_parsed;
    }

    size_t n_skipped() const {
        return _n_skipped;
    }

    bool is_complete() const {
        return _is_complete;
    }
}; // class FastxParser


template <class ParserType = FastxParser<>>
class SplitPairedReader {

    typedef ParserType                     parser_type;
    typedef typename parser_type::alphabet alphabet;

    std::shared_ptr<parser_type> left_parser;
    std::shared_ptr<parser_type> right_parser;
    bool                         _force_name_match;
    bool                         _strict;
    uint64_t                     _n_skipped;

public:

    typedef RecordPair value_type;

    SplitPairedReader(const std::string &left,
                      const std::string &right,
                      bool strict = false,
                      uint32_t min_length = 0,
                      bool force_name_match = false)
        :  _force_name_match(force_name_match),
          _strict(strict),
          _n_skipped(0) {
        
        left_parser = parser_type::build(left, strict, min_length);
        right_parser = parser_type::build(right, strict, min_length);
    }

    static std::shared_ptr<SplitPairedReader<ParserType>> build(const std::string &left,
                                                                const std::string &right,
                                                                bool strict = false,
                                                                uint32_t min_length = 0,
                                                                bool force_name_match = false) {
        return std::make_shared<SplitPairedReader<ParserType>>(left, right, strict, min_length, force_name_match);
    }

    bool is_complete() const {
        if (left_parser->is_complete() != right_parser->is_complete()) {
            throw InvalidPairedStream();
        }
        return left_parser->is_complete();
    }

    RecordPair next() {

        if (is_complete()) {
            throw std::out_of_range("No more reads available.");
        }

        std::optional<Record> left, right;
        std::exception_ptr left_exc_ptr, right_exc_ptr;
        
        try {
            left = this->left_parser->next();
        } catch (InvalidSequence) {
            left_exc_ptr = std::current_exception();
        } catch (InvalidRecord) {
            left_exc_ptr = std::current_exception();
        }

        try {
            right = this->right_parser->next();
        } catch (InvalidSequence) {
            right_exc_ptr = std::current_exception();
        } catch (InvalidRecord) {
            right_exc_ptr = std::current_exception();
        }

        if (_strict) {
            if (left_exc_ptr) {
                std::rethrow_exception(left_exc_ptr);
            }
            if (right_exc_ptr) {
                std::rethrow_exception(right_exc_ptr);
            }
        }

        if (this->_force_name_match && left && right) {
            if (!check_is_pair(left.value().name, right.value().name)) {
                if (_strict) {
                    throw InvalidRecordPair(this->left_parser->n_parsed(),
                                            left.value().name,
                                            right.value().name);
                } else {
                    _n_skipped += 2;
                    return {{}, {}};
                }
            }
        }

        return std::make_pair(left, right);
    }

    uint64_t n_skipped() const {
        return _n_skipped + left_parser->n_skipped() + right_parser->n_skipped();
    }

    uint64_t n_parsed() const {
        return left_parser->n_parsed() + right_parser->n_parsed();
    }
};


extern template class goetia::FastxParser<goetia::DNA_SIMPLE>;
extern template class goetia::FastxParser<goetia::DNAN_SIMPLE>;
extern template class goetia::FastxParser<goetia::IUPAC_NUCL>;

extern template class goetia::SplitPairedReader<goetia::FastxParser<goetia::DNA_SIMPLE>>;
extern template class goetia::SplitPairedReader<goetia::FastxParser<goetia::DNAN_SIMPLE>>;
extern template class goetia::SplitPairedReader<goetia::FastxParser<goetia::IUPAC_NUCL>>;

} // namespace goetia

#endif 
