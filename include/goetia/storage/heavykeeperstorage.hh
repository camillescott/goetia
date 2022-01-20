/**
 * (c) Camille Scott, 2021
 * File   : heavykeeperstorage.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 08.11.2021
 */

#ifndef GOETIA_HKSTORAGE_HH
#define GOETIA_HKSTORAGE_HH

#include "goetia/sketches/sketch/include/sketch/hk.h"
#include "goetia/storage/storage.hh"

namespace goetia::storage {

class HKStorage;

template<>
struct StorageTraits<HKStorage> {
    static constexpr bool is_probabilistic = true;
    static constexpr bool is_counting = true;
    typedef std::tuple<size_t> params_type;
};

}

#endif
