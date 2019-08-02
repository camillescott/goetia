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

#include "boink/assembly.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/dbg.hh"
#include "boink/cdbg/udbg.hh"
#include "boink/cdbg/compactor.hh"

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

template <class GraphType>
struct UStreamingCompactor {

    using ShifterType     = typename GraphType::shifter_type;
    using TraversalType   = Traverse<GraphType>;
    using TraverserType   = typename TraversalType::dBG;
    using MinimizerType   = WKMinimizer<ShifterType>;
    using UnikmerIterType = typename hashing::KmerIterator<hashing::UKHShifter>;


    using cDBGType      = uDBG<GraphType>;
    using UnitigTip     = typename cDBGType::UnitigTip;
    using UnitigTag     = typename cDBGType::Tag;
    using State         = typename TraversalType::State;

    typedef TipType::Type MarkerType;
    struct TipMarker {
        kmer_t              kmer;
        std::vector<kmer_t> neighbors;
        direction_t         side;
    };

    struct Segment {
        kmer_t                                ltip, rtip;
        std::vector<kmer_t>                   lneighbors, rneighbors;
        std::vector<hashing::PartitionedHash> hashes;
    };


    class Compactor : public TraversalType::dBG,
                      public events::EventNotifier {

    public:
 
        using TraversalType::dBG::filter_nodes;
        using TraversalType::dBG::find_left_kmers;
        using TraversalType::dBG::find_right_kmers;
        using TraversalType::dBG::gather_left;
        using TraversalType::dBG::gather_right;
        using TraversalType::dBG::traverse_left;
        using TraversalType::dBG::traverse_right;
        using TraversalType::dBG::get_decision_neighbors;

        typedef ShifterType   shifter_type;
        typedef GraphType     graph_type;
        typedef MinimizerType minimizer_type;
        typedef cDBGType      cdbg_type;


        // Graphs: underling dBG, cDBG tip graph
        shared_ptr<GraphType>                dbg;
        shared_ptr<typename cDBGType::Graph> cdbg;

        // Universal k-mer hitting set
        shared_ptr<hashing::UKHS>            ukhs;
        const uint16_t                       unikmer_k;
        hashing::UKHShifter                  unikmer_partitioner;

        Compactor(std::shared_ptr<GraphType> dbg,
                  std::shared_ptr<hashing::UKHS> ukhs)
            : TraversalType::dBG(dbg->K()),
              EventNotifier(),
              dbg(dbg),
              ukhs(ukhs),
              unikmer_k(ukhs->K()),
              unikmer_partitioner(dbg->K(), ukhs->K(), ukhs)
        {
            this->cdbg = std::make_shared<typename cDBGType::Graph>(dbg);
        }

        static std::shared_ptr<Compactor> build(std::shared_ptr<GraphType> dbg,
                                                std::shared_ptr<hashing::UKHS> ukhs) {
            return std::make_shared<Compactor>(dbg, ukhs);
        }

        void find_segment_seeds(const std::string&                     sequence,
                                std::deque<TraverserType>&             seeds,
                                std::deque<std::pair<size_t, size_t>>& positions,
                                std::set<hash_t>&                      new_kmers,
                                std::deque<hashing::PartitionedHash>&  ordered_pkmers) {
            
            UnikmerIterType kmers(sequence, &unikmer_partitioner);

            hashing::PartitionedHash prev_hash, cur_hash;
            size_t pos = 0, start_pos = ULLONG_MAX;
            bool cur_new = false, prev_new = false, cur_seen = false, prev_seen = false;
            
            while(!kmers.done()) {
                cur_hash = kmers.next();
                cur_new = this->dbg->query(cur_hash.first) == 0;
                cur_seen = new_kmers.count(cur_hash.first);

                if(cur_new) {
                    new_kmers.insert(cur_hash.first);
                    ordered_pkmers.push_back(cur_hash);
                    if (!cur_seen) {
                        if(!prev_new || prev_seen) {
                            pdebug("old -> new (pos=" << pos << ")");
                            start_pos = pos;
                            seeds.emplace_back(typename TraversalType::dBG(this->_K));
                            seeds.back().set_cursor(sequence.c_str() + pos);
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
                                           std::deque<TraverserType>&             seeds,
                                           std::deque<std::pair<size_t, size_t>>& positions,
                                           std::set<hash_t>&                      new_kmers,
                                           std::deque<hashing::PartitionedHash>&  ordered_pkmers) {

            std::deque<TipMarker> markers;

            for (TraverserType& seed : seeds) {
                
                auto coords = positions.pop_front();
                size_t pos = coords.first;
                size_t shift_pos = pos + this->_K - 1;

                std::deque<hashing::PartitionedHash> segment_pkmers;
                std::vector<kmer_t> left_neighbors, right_neighbors;
                kmer_t              current;

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
