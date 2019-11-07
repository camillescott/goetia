/**
 * (c) Camille Scott, 2019
 * File   : meta.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 05.11.2019
 */

#ifndef BOINK_METAUTIL_HH
#define BOINK_METAUTIL_HH

#include <string>

namespace boink {

struct Versioned {
    static const unsigned int VERSION;
};

template <typename T>
struct Tagged : public Versioned {
    static const std::string NAME;
};

}

#endif
