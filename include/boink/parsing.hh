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
#include "utils/stringutils.h"

#include "boink/boink.hh"

using namespace oxli;
using namespace oxli:: read_parsers;

using namespace utils;

namespace boink {


struct ReadBundle {
    
    bool has_left;
    bool has_right;

    Read left;
    Read right;

};


bool check_char(const char c, const std::string against) {
    for (auto ag : against) {
        if (c == ag) {
            return true;
        }
    }
    return false;
}


std::pair<std::string, std::string> split_on_first(const std::string& name,
                                    const std::string delims=" \t") {
    std::string left = "";
    std::string right = "";
    for (size_t i = 0; i < name.length(); ++i) {
        const char c = name[i];
        if (left == "") {
            if (check_char(c, delims)) {
                left = name.substr(0, i);
            }
        } else {
            if (!check_char(c, delims)) {
                right = name.substr(i, name.length());
                break;
            }
        }
    }
    if (left == "") {
        left = name;
    }

    return std::make_pair(left, right);
}


bool check_is_pair(const std::string& left, const std::string& right) {
    auto left_split = split_on_first(left);
    auto right_split = split_on_first(right);
    
    if (StringUtils::endsWith(left_split.first, "/1") &&
        StringUtils::endsWith(right_split.first, "/2")) {
        
        auto left_sub = split_on_first(left_split.first, "/");
        auto right_sub = split_on_first(right_split.first, "/");

        if (left_sub.first != ""  && left_sub.first == right_sub.first) {
            return true;
        }
    } else if (left_split.first == right_split.first &&
               StringUtils::endsWith(left_split.second, "1:") &&
               StringUtils::endsWith(right_split.second, "2:")) {
        return true;
    } else if (left_split.first == right_split.first &&
               StringUtils::endsWith(left_split.second, "/1") &&
               StringUtils::endsWith(right_split.second, "/2")) {
        auto left_sub = split_on_first(left_split.second, "/");
        auto right_sub = split_on_first(right_split.second, "/");

        if (left_sub.first != "" && left_sub.first == right_sub.first) {
            return true;
        }
    }
    return false;
}


bool check_is_left(const std::string& name) {
    auto split = split_on_first(name);
    if (StringUtils::endsWith(split.first, "/1") ||
        StringUtils::endsWith(split.second, "1:") ||
        StringUtils::endsWith(split.second, "/1")) {
        return true;
    }
    return false;
}


bool check_is_right(const std::string& name) {
    auto split = split_on_first(name);
    if (StringUtils::endsWith(split.first, "/2") ||
        StringUtils::endsWith(split.second, "2:") ||
        StringUtils::endsWith(split.second, "/2")) {
        return true;
    }
    return false;
}

void filter_length(ReadBundle& bundle, uint32_t length) {
    if (bundle.left.sequence.length() < length) {
        bundle.has_left = false;
    }
    if (bundle.right.sequence.length() < length) {
        bundle.has_right = false;
    }
}


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

}

#endif
