/* report_types.hh 
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_REPORT_TYPES_HH
#define BOINK_REPORT_TYPES_HH

#include <cstdint>

namespace boink {
namespace reporting {

struct StreamingCompactorReport {
    uint64_t n_full;
    uint64_t n_tips;
    uint64_t n_islands;
    uint64_t n_trivial;
    uint64_t n_circular;
    uint64_t n_loops;
    uint64_t n_dnodes;
    uint64_t n_unodes;
    uint64_t n_updates;
    uint64_t n_splits;
    uint64_t n_merges;
    uint64_t n_extends;
    uint64_t n_clips;
    uint64_t n_deletes;
    uint64_t n_circular_merges;
    uint64_t n_tags;
    uint64_t n_unique;
    double   estimated_fp;
};


}
}

#endif
