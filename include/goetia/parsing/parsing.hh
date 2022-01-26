/* parsing.hh -- parsers
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef GOETIA_PARSING_HH
#define GOETIA_PARSING_HH


#include <algorithm>
#include <iostream>
#include <optional>
#include <string>
#include <utility>

#include "goetia/goetia.hh"


namespace goetia {


struct Record {
    std::string name;
    std::string sequence;
    std::string quality;

    inline void reset()
    {
        name.clear();
        sequence.clear();
        quality.clear();
    }

    inline void write_fastx(std::ostream& output) const
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

    friend inline std::ostream& operator<<(std::ostream& o, const Record& record) {
        o << "<Sequence name=" << record.name
          << " seq=" << record.sequence
          << ">";
        return o;
    }
};


typedef std::pair<std::optional<Record>, std::optional<Record>> RecordPair;

bool check_char(const char c, const std::string against);

std::pair<std::string, std::string> split_on_first(const std::string& name,
                                                   const std::string delims);

bool check_is_pair(const std::string& left, const std::string& right);

bool check_is_left(const std::string& name);

bool check_is_right(const std::string& name);

}

#endif
