/* goetia.hh
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "goetia/parsing/parsing.hh"
#include "goetia/utils/stringutils.h"

#include <utility>
#include <string>


using namespace utils;

namespace goetia {
namespace parsing {

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


bool check_is_pair(const std::string& left,
                   const std::string& right) {

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

}
}
