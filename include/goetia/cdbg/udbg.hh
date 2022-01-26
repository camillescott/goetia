/**
 * (c) Camille Scott, 2019
 * File   : udbg.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 21.06.2019
 */


#ifndef GOETIA_UDBG_HH
#define GOETIA_UDBG_HH

// save diagnostic state
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wchar-subscripts"
#include "goetia/storage/sparsepp/spp.h"
#pragma GCC diagnostic pop

#include "goetia/goetia.hh"
#include "goetia/meta.hh"

#include "goetia/dbg.hh"
#include "goetia/traversal/unitig_walker.hh"
#include "goetia/hashing/kmeriterator.hh"
#include "goetia/hashing/rollinghashshifter.hh"
#include "goetia/hashing/ukhs.hh"
#include "goetia/storage/storage_types.hh"

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

namespace goetia {


template <class StorageType>
struct uDBG {

    typedef dBG<StorageType, CanUnikmerShifter> graph_type;

    // inject dependent typename boilerplate: see goetia/meta.hh
    _goetia_model_typedefs_from_graphtype(graph_type);
    _goetia_walker_typedefs_from_graphtype(graph_type);


    struct UnitigTip {
        std::string             kmer;
        std::vector<value_type> neighbors;
        bool                    position;
        value_type              partner;

        UnitigTip() {}

        UnitigTip(std::string kmer,
                  bool position)
            : kmer(kmer),
              position(position) {

        }
    };

    struct Tag {
        value_type  left;
        value_type  right;
        uint64_t    partition;
    };

    class Graph {

        typedef spp::sparse_hash_map<value_type, UnitigTip*> tip_map_t;
        typedef typename tip_map_t::const_iterator           tip_map_iter_t;

        typedef spp::sparse_hash_map<value_type, Tag>        tag_map_t;
        typedef typename tag_map_t::const_iterator           tag_map_iter_t;

    protected:

        std::shared_ptr<graph_type> dbg;

        tip_map_t  left_tips;
        tip_map_t  right_tips;
        tag_map_t  tags;

        uint64_t   _n_decision_nodes;
        uint64_t   _n_updates;


    public:

        const uint16_t K;

        Graph(std::shared_ptr<graph_type> dbg)
            : K(dbg->K),
              dbg(dbg) {
        }

        UnitigTip * query_left_tips(const hash_type& tip_hash) {
            auto search = left_tips.find(tip_hash.value());
            if (search != left_tips.end()) {
                return search->second;
            }
            return nullptr;
        }

        UnitigTip * query_right_tips(const hash_type& tip_hash) {
            auto search = right_tips.find(tip_hash.value());
            if (search != right_tips.end()) {
                return search->second;
            }
            return nullptr;
        }

        state_type left_traverse_to_tip(const hash_type& start,
                                        UnitigTip *      result) {

            result = nullptr;
            auto search = tags.find(start.value());
            if (search == tags.end()) {
                return state_type::BAD_SEED;
            }
            auto cur = search->second;

            std::set<value_type> seen;
            seen.insert(start.value());
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

extern template class uDBG<BitStorage>;
extern template class uDBG<ByteStorage>;
extern template class uDBG<NibbleStorage>;
extern template class uDBG<QFStorage>;
extern template class uDBG<SparseppSetStorage>;



}

#undef pdebug
#endif
