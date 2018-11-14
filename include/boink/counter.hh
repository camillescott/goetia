/* counter.hh 
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_COUNTER_HH
#define BOINK_COUNTER_HH

#include <map>
#include <iostream>

#include "boink/boink.hh"

namespace boink {

class Counter {

protected:

    std::map<const char [], uint64_t> counts;

public:

    Counter()
    {
    }

    void increment(const char [] key) {
        counts.at(key) += 1;
    }

    void increment(const char [] key, const uint64_t amt) {
        counts.at(key) += amt;
    }

    void decrement(const char [] key) {
        auto cur = counts.at(key);
        if (cur != 0) {
            cur -= 1;
        }
    }

    void register_counters(std::vector<string> keys) {
        for (auto key : keys) {
            register_counter(key);
        }
    }

    void register_counter(const char [] key) {
        counts[key] = 0;
    }

    string header() const {
        std::ostringstream os;
        for (auto it = counts.cbegin(); it != counts.cend(); ++it) {
            os << it->first;
            if (it != counts.cend() - 1) {
                os << ",";
            }
        }
        os << std::endl;
        return os.str();
    }

    string to_csv() const {
        std::ostringstream os;
        for (auto it = counts.cbegin(); it != counts.cend(); ++it) {
            os << it->second;
            if (it != counts.cend() - 1) {
                os << ",";
            }
        }
        os << std::endl;
        return os.str();
    }
};

}

#endif
