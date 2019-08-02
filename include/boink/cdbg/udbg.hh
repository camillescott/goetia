/**
 * (c) Camille Scott, 2019
 * File   : udbg.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 21.06.2019
 */
/* udbg.hh -- compact de Bruijn Graph
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_UDBG_HH
#define BOINK_UDBG_HH

// save diagnostic state
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wchar-subscripts"
#include <gfakluge.hpp>
#include "boink/storage/sparsepp/spp.h"
#pragma GCC diagnostic pop

#include "boink/boink.hh"

#include "boink/events.hh"
#include "boink/assembly.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/kmers/kmerclient.hh"
#include "boink/minimizers.hh"
#include "boink/storage/storage.hh"

# ifdef DEBUG_CDBG
#   define pdebug(x) do { std::ostringstream stream; \
                          stream << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          std::cerr << stream.str(); \
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif

namespace boink {
namespace cdbg {

using boink::hashing::hash_t;
using boink::hashing::shift_t;
using boink::hashing::kmer_t;
using boink::hashing::KmerIterator;


template <class GraphType>
struct uDBG {

    protected:

        using ShifterType   = typename GraphType::shifter_type;
        using TraversalType = Traverse<GraphType>;
        using MinimizerType = typename WKMinimizer<ShifterType>::Minimizer;

    public:

    struct Junction {
        hash_t src;
        hash_t dst;
    };

    struct UnitigTip {
        std::string         kmer;
        std::vector<hash_t> neighbors;
        direction_t         position;
        hash_t              partner;

        UnitigTip() {}

        UnitigTip(std::string kmer,
                  direction_t position)
            : kmer(kmer),
              position(position) {

        }
    };

    struct Tag {
        hash_t   left;
        hash_t   right;
        uint64_t unikmer;
    };

    typedef GraphType             graph_type;
    typedef ShifterType           shifter_type;
    typedef TraversalType         traverser_type;
    typedef MinimizerType         minimizer_type;
    typedef TraversalState::State state_type;

    class Graph : public kmers::KmerClient,
                  public events::EventNotifier {

        typedef spp::sparse_hash_map<hash_t, UnitigTip*> tip_map_t;
        typedef typename tip_map_t::const_iterator       tip_map_iter_t;

        typedef spp::sparse_hash_map<hash_t, Tag>        tag_map_t;
        typedef typename tag_map_t::const_iterator       tag_map_iter_t;

    protected:

        std::shared_ptr<GraphType> dbg;

        tip_map_t  left_tips;
        tip_map_t  right_tips;
        tag_map_t  tags;

        uint64_t   _n_decision_nodes;
        uint64_t   _n_updates;


    public:

        Graph(std::shared_ptr<GraphType> dbg)
            : KmerClient(dbg->K()),
              dbg(dbg) {
        }

        UnitigTip * query_left_tips(hash_t tip_hash) {
            auto search = left_tips.find(tip_hash);
            if (search != left_tips.end()) {
                return search->second;
            }
            return nullptr;
        }

        UnitigTip * query_right_tips(hash_t tip_hash) {
            auto search = right_tips.find(tip_hash);
            if (search != right_tips.end()) {
                return search->second;
            }
            return nullptr;
        }

        state_type left_traverse_to_tip(hash_t      start,
                                        UnitigTip * result) {

            result = nullptr;
            auto search = tags.find(start);
            if (search == tags.end()) {
                return state_type::BAD_SEED;
            }
            auto cur = search->second;

            std::set<hash_t> seen;
            seen.insert(start);
            while(true) {
                search = tags.find(cur.left);
                if (search == tags.end()) {
                    auto tip_search = left_tips.find(cur.left);
                    if (tip_search != left_tips.end()) {
                        result = tip_search->second;
                        return state_type::STOP_FWD;
                    } else {
                        return state_type::GRAPH_ERROR;
                    }
                }
                if (seen.count(cur.left)) {
                    return state_type::STOP_SEEN;
                }
                seen.insert(cur.left);
                cur = search->second;
            }
        }

    };

};

}
}

#undef pdebug
#endif
