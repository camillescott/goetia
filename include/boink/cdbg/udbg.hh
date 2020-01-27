/**
 * (c) Camille Scott, 2019
 * File   : udbg.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 21.06.2019
 */


#ifndef BOINK_UDBG_HH
#define BOINK_UDBG_HH

// save diagnostic state
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wchar-subscripts"
#include "boink/storage/sparsepp/spp.h"
#pragma GCC diagnostic pop

#include "boink/boink.hh"
#include "boink/meta.hh"

#include "boink/events.hh"
#include "boink/dbg.hh"
#include "boink/traversal.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/storage/storage_types.hh"

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


template <class StorageType>
struct uDBG {

    typedef dBG<StorageType, hashing::CanUnikmerShifter> graph_type;

    // inject dependent typename boilerplate: see boink/meta.hh
    _boink_model_typedefs_from_graphtype(graph_type);
    _boink_walker_typedefs_from_graphtype(graph_type);


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

    class Graph : public events::EventNotifier {

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

extern template class cdbg::uDBG<storage::BitStorage>;
extern template class cdbg::uDBG<storage::ByteStorage>;
extern template class cdbg::uDBG<storage::NibbleStorage>;
extern template class cdbg::uDBG<storage::QFStorage>;
extern template class cdbg::uDBG<storage::SparseppSetStorage>;



}
}

#undef pdebug
#endif
