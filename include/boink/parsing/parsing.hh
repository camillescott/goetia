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


#include <algorithm>
#include <string>
#include <utility>

#include "boink/boink.hh"


namespace boink {
namespace parsing {

unsigned char _to_valid_dna(const unsigned char c);

struct Read {
    std::string name;
    std::string description;
    std::string sequence;
    std::string quality;
    std::string cleaned_seq;

    inline void reset()
    {
        name.clear();
        description.clear();
        sequence.clear();
        quality.clear();
        cleaned_seq.clear();
    }

    inline void write_fastx(std::ostream& output)
    {
        if (quality.length() != 0) {
            output << "@" << name << '\n'
                   << sequence << '\n'
                   << "+" << '\n'
                   << quality << '\n';
        } else {
            output << ">" << name << '\n'
                   << sequence << '\n';
        }
    }

    // Compute cleaned_seq from sequence. Call this after changing sequence.
    inline void set_clean_seq()
    {
        cleaned_seq = std::string(sequence.size(), 0);
        std::transform(sequence.begin(), sequence.end(), cleaned_seq.begin(),
                       _to_valid_dna);
    }
};


typedef std::pair<Read, Read> ReadPair;


struct ReadBundle {
    
    bool has_left;
    bool has_right;

    Read left;
    Read right;

};

bool check_char(const char c, const std::string against);

std::pair<std::string, std::string> split_on_first(const std::string& name,
                                                   const std::string delims);

bool check_is_pair(const std::string& left, const std::string& right);

bool check_is_left(const std::string& name);

bool check_is_right(const std::string& name);

void filter_length(ReadBundle& bundle, uint32_t length);

}
}

#endif
