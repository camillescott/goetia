/* parsing.hh -- parsers
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_PARSING_HH
#define BOINK_PARSING_HH

#include <string>
#include <utility>

#include "oxli/read_parsers.hh"

#include "boink/boink.hh"


namespace boink {


struct ReadBundle {
    
    bool has_left;
    bool has_right;

    oxli::read_parsers::Read left;
    oxli::read_parsers::Read right;

};


bool check_char(const char c, const std::string against);

std::pair<std::string, std::string> split_on_first(const std::string& name,
                                                   const std::string delims);

bool check_is_pair(const std::string& left, const std::string& right);

bool check_is_left(const std::string& name);

bool check_is_right(const std::string& name);

void filter_length(ReadBundle& bundle, uint32_t length);

template <class ParserType = oxli::read_parsers::FastxReader>
class SplitPairedReader {

    oxli::read_parsers::ReadParserPtr<ParserType> left_parser;
    oxli::read_parsers::ReadParserPtr<ParserType> right_parser;
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
        
        left_parser = oxli::read_parsers::get_parser<ParserType>(left);
        right_parser = oxli::read_parsers::get_parser<ParserType>(right);
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
        } catch (oxli::read_parsers::NoMoreReadsAvailable) {
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

}

#endif
