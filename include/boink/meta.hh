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
#include <utility>
#include <vector>

#define _boink_model_typedefs_from_shiftertype(shifter_type) \
    typedef hashing::HashExtender<shifter_type> extender_type; \
    typedef typename shifter_type::alphabet     alphabet; \
    typedef typename shifter_type::hash_type    hash_type; \
	typedef typename hash_type::value_type      value_type; \
    typedef typename shifter_type::kmer_type    kmer_type; \
    template<hashing::Direction_t D> using shift_type = hashing::ShiftModel<hash_type, D>; \
    typedef std::pair<std::vector<kmer_type>, \
                      std::vector<kmer_type>>                      neighbor_pair_type;\
    typedef std::pair<std::vector<shift_type<hashing::DIR_LEFT>>,\
                      std::vector<shift_type<hashing::DIR_RIGHT>>> shift_pair_type;

#define _boink_model_typedefs_from_graphtype(graph_type) \
    typedef typename graph_type::shifter_type   shifter_type; \
    _boink_model_typedefs_from_shiftertype(shifter_type) \

    
#define _boink_walker_typedefs_from_graphtype(graph_type) \
    typedef dBGWalker<graph_type> walker_type; \
    template<hashing::Direction_t D>  using walk_type = typename walker_type::template Walk<D>;\
    typedef typename walker_type::walk_pair_type walk_pair_type; \
    typedef typename walker_type::State state_type;

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
