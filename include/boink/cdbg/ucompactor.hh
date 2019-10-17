/**
 * (c) Camille Scott, 2019
 * File   : ucompactor.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 27.06.2019
 */
/* compactor.hh -- streaming dBG compactor
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_UCOMPACTOR_HH
#define BOINK_UCOMPACTOR_HH

#include <assert.h>
#include <climits>
#include <cstdint>

#include "boink/traversal.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/dbg.hh"
#include "boink/cdbg/udbg.hh"

# ifdef DEBUG_CPTR
#   define pdebug(x) do { std::cerr << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif

namespace boink {
namespace cdbg {

namespace TipType {
    enum Type {
        LEFT,
        RIGHT,
        LEFT_FLANK,
        RIGHT_FLANK
    };
}

template <class StorageType>
struct UStreamingCompactor {

    using dbg_type          = dBG<StorageType,
                                  hashing::UKHS::LazyShifter>;
    using storage_type      = StorageType;
    using shifter_type      = typename dbg_type::shifter_type;
    using traversal_type    = Traverse<dbg_type>;
    using traverser_type    = typename traversal_type::dBG;
    using kmer_iter_type    = typename hashing::KmerIterator<hashing::UKHS::LazyShifter>;
    using hash_type         = typename dbg_type::hash_type;
    using kmer_type         = typename dbg_type::kmer_type;
    using shift_type        = typename dbg_type::shift_type;

    using cdbg_type     = uDBG<storage_type>;
    using UnitigTip     = typename cdbg_type::UnitigTip;
    using UnitigTag     = typename cdbg_type::Tag;
    using State         = typename traversal_type::State;

    typedef TipType::Type MarkerType;
    struct TipMarker {
        kmer_type              kmer;
        std::vector<kmer_type> neighbors;
        direction_t            side;
    };

    struct Segment {
        kmer_type              ltip, rtip;
        std::vector<kmer_type> lneighbors, rneighbors;
        std::vector<hash_type> hashes;
    };


    class Compactor : public traverser_type,
                      public events::EventNotifier {

    public:
 
        using traverser_type::filter_nodes;
        using traverser_type::find_left_kmers;
        using traverser_type::find_right_kmers;
        using traverser_type::gather_left;
        using traverser_type::gather_right;
        using traverser_type::traverse_left;
        using traverser_type::traverse_right;
        using traverser_type::get_decision_neighbors;



        // Graphs: underling dBG, cDBG tip graph
        shared_ptr<dbg_type>                 dbg;
        shared_ptr<typename cdbg_type::Graph> cdbg;

        // Universal k-mer hitting set
        shared_ptr<hashing::UKHS::Map>       ukhs;
        const uint16_t                       unikmer_k;

        Compactor(std::shared_ptr<dbg_type> dbg,
                  std::shared_ptr<hashing::UKHS::Map> ukhs)
            : traverser_type(dbg->K(), ukhs->K(), ukhs),
              EventNotifier(),
              dbg(dbg),
              ukhs(ukhs),
              unikmer_k(ukhs->K())
        {
            this->cdbg = std::make_shared<typename cdbg_type::Graph>(dbg);
        }

        static std::shared_ptr<Compactor> build(std::shared_ptr<dbg_type> dbg,
                                                std::shared_ptr<hashing::UKHS::Map> ukhs) {
            return std::make_shared<Compactor>(dbg, ukhs);
        }

        void find_segment_seeds(const std::string&                     sequence,
                                //std::deque<traverser_type>&            seeds,
                                std::deque<std::pair<size_t, size_t>>& positions,
                                std::set<hash_type>&                   new_kmers,
                                std::deque<hash_type>&                 ordered_pkmers) {
            
            kmer_iter_type kmers(sequence, this);

            hash_type prev_hash, cur_hash;
            size_t pos = 0, start_pos = ULLONG_MAX;
            bool cur_new = false, prev_new = false, cur_seen = false, prev_seen = false;
            
            while(!kmers.done()) {
                cur_hash = kmers.next();
                cur_new = this->dbg->query(cur_hash.hash) == 0;
                cur_seen = new_kmers.count(cur_hash.hash);

                if(cur_new) {
                    new_kmers.insert(cur_hash.hash);
                    ordered_pkmers.push_back(cur_hash);
                    if (!cur_seen) {
                        if(!prev_new || prev_seen) {
                            pdebug("old -> new (pos=" << pos << ")");
                            start_pos = pos;
                            //seeds.emplace_back({this->_K, this->unikmer_k, this->ukhs});
                            //seeds.back().set_cursor(sequence.c_str() + pos);
                        }
                    }
                } else if (start_pos != ULLONG_MAX && (!cur_new || cur_seen)) {
                    pdebug("new -> old, was_seen=" << cur_seen);
                    positions.emplace_back(start_pos, pos - 1);
                    start_pos = ULLONG_MAX;
                } 

                ++pos;
                prev_hash = cur_hash;
                prev_new = cur_new;
                prev_seen = cur_seen;

            }

            if (cur_new && !cur_seen) {
                pdebug("sequence ended on new k-mer");
                positions.emplace_back(start_pos, pos - 1);
            }
        }
        
        /*
        std::deque<Segment> build_segments(const std::string&                     sequence,
                                           std::deque<traverser_type>&             seeds,
                                           std::deque<std::pair<size_t, size_t>>& positions,
                                           std::set<hash_t>&                      new_kmers,
                                           std::deque<hash_type>&  ordered_pkmers) {

            std::deque<TipMarker> markers;

            for (traverser_type& seed : seeds) {
                
                auto coords = positions.pop_front();
                size_t pos = coords.first;
                size_t shift_pos = pos + this->_K - 1;

                std::deque<hash_type> segment_pkmers;
                std::vector<kmer_type> left_neighbors, right_neighbors;
                kmer_type              current;

                // beginning of seed: find left anchor
                left_neighbors = seed->find_left_kmers(dbg.get(),
                                                       new_kmers);
                right_neighbors = seed->find_right_kmers(dbg.get(),
                                                         new_kmers);
                seed->get_cursor(current);
                
                // search for right anchor or end of seed segment           
                while (true) {

                    right_neighbors = seed->find_right_kmers(dbg.get(),
                                                             new_kmers);
                    if (right_neighbors.size() > 1) {
                        seed->get_cursor(right_end);
                        
                    }
                }


            }
        }
        */

        
                         
    };


};

}
}

#undef pdebug
#endif
