/* compactor.hh -- streaming dBG compactor
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_SOLID_COMPACTOR_HH
#define BOINK_SOLID_COMPACTOR_HH

#include <assert.h>

#include "boink/assembly.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/dbg.hh"
#include "boink/cdbg/cdbg.hh"
#include "boink/minimizers.hh"
#include "boink/event_types.hh"
#include "boink/reporting/report_types.hh"
#include "boink/storage/bytestorage.hh"

#include "boink/cdbg/compactor.hh"


namespace boink {
namespace cdbg {

using namespace boink::hashing;
using namespace boink::storage;


# ifdef DEBUG_CPTR
#   define pdebug(x) do { std::cerr << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif


template<class GraphType>
class SolidStreamingCompactor : public events::EventNotifier {

private:

    std::unique_ptr<dBG<ByteStorage,
                        typename GraphType::shifter_type>> abund_filter;

public:

    std::shared_ptr<StreamingCompactor<GraphType>> compactor;
    std::shared_ptr<GraphType>                     dbg;

    unsigned int                                   min_abund;

    SolidStreamingCompactor(shared_ptr<StreamingCompactor<GraphType>> compactor,
                            unsigned int                              min_abund,
                            uint64_t                                  abund_table_size,
                            uint16_t                                  n_abund_tables)
        : EventNotifier (),
          compactor     (compactor),
          dbg           (compactor->dbg),
          min_abund     (min_abund)
    {
        abund_filter = std::make_unique<dBG<ByteStorage,
                                            typename GraphType::shifter_type>>(dbg->K(),
                                                                               abund_table_size,
                                                                               n_abund_tables);
    }

    std::vector<std::pair<size_t, size_t>> find_solid_segments(const std::string& sequence) {
        std::vector<hash_t>  hashes;
        std::vector<count_t> counts;
        abund_filter->insert_sequence(sequence, hashes, counts);

        std::vector<std::pair<size_t, size_t>> segments;
        auto hashes_iter = hashes.begin();
        auto counts_iter = counts.begin();
        size_t start, end, pos = 0;
        hash_t cur_hash;
        bool   prev_solid = false, cur_solid = false;
        while(hashes_iter != hashes.end()) {
            cur_hash  = *hashes_iter;
            cur_solid = (*counts_iter) >= min_abund;
            if (cur_solid && !prev_solid) {
                start = pos;
            } else if (prev_solid && !cur_solid) {
                end = pos + compactor->K() - 1;
                segments.emplace_back(start, end);
            }
            
            ++pos;
            ++hashes_iter;
            ++counts_iter;
            prev_solid = cur_solid;
        }

        if (cur_solid) {
            end = pos + compactor->K() - 1;
            segments.emplace_back(start, end);
        }

        return segments;
    }

    void update_sequence(const std::string& sequence) {

        auto solid_segments = find_solid_segments(sequence);
        for (auto segment : solid_segments) {
            auto length = segment.second - segment.first;
            compactor->update_sequence(sequence.substr(segment.first,
                                                       length));
        }
    }
          
                                
};

}
}

#endif
