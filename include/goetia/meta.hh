/**
 * (c) Camille Scott, 2019
 * File   : meta.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 05.11.2019
 */

#ifndef GOETIA_METAUTIL_HH
#define GOETIA_METAUTIL_HH

#include <string>
#include <utility>
#include <vector>

#include "goetia/nameof.hpp"

#define _goetia_model_typedefs_from_shiftertype(shifter_type) \
    typedef typename hashing::extender_selector<shifter_type>::type extender_type; \
    typedef typename shifter_type::alphabet          alphabet; \
    typedef typename shifter_type::hash_type         hash_type; \
	typedef typename hash_type::value_type           value_type; \
    typedef typename shifter_type::kmer_type         kmer_type; \
    template<bool Dir> using shift_type = hashing::Shift<hash_type, Dir>; \
    typedef std::pair<std::vector<kmer_type>, \
                      std::vector<kmer_type>>         neighbor_pair_type;\
    typedef std::pair<std::vector<shift_type<hashing::DIR_LEFT>>,\
                      std::vector<shift_type<hashing::DIR_RIGHT>>> shift_pair_type;

#define _goetia_model_typedefs_from_graphtype(graph_type) \
    typedef typename graph_type::shifter_type   shifter_type; \
    _goetia_model_typedefs_from_shiftertype(shifter_type) \

    
#define _goetia_walker_typedefs_from_graphtype(walker_type) \
    template<bool Dir> \
        using walk_type = typename walker_type::template Walk<Dir>;\
    typedef typename walker_type::walk_pair_type walk_pair_type; \
    typedef typename walker_type::State state_type;

namespace goetia {

static constexpr size_t GLOBAL_ABI_VERSION = 0;

template <size_t V = 0>
struct Versioned {
    //static constexpr size_t GLOBAL_ABI_VERSION = GLOBAL_ABI_VERSION;
    static constexpr size_t OBJECT_ABI_VERSION = V;

    static const std::string version_string() {
        return std::to_string(OBJECT_ABI_VERSION);
    }

    static const char * version_binary() {
        return reinterpret_cast<const char *>(&OBJECT_ABI_VERSION);
    }
};

template <typename T,
          size_t V = 0>
struct Tagged : public Versioned<V> {
    using Versioned<V>::OBJECT_ABI_VERSION;
    using Versioned<V>::version_string;
    using Versioned<V>::version_binary;

    static constexpr auto NAME = NAMEOF_TYPE(T);

    static const std::string name_string() {
        return std::string(NAMEOF_TYPE(T));
    }
};


}

#endif
