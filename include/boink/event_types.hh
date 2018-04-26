/* report_types.hh 
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef EVENT_TYPES_HH
#define EVENT_TYPES_HH

namespace boink {
namespace event_types {

enum event_t {
    // EventListener control events
    MSG_EXIT_THREAD,
    MSG_TIMER,

    // Data events
    MSG_WRITE_CDBG_STATS,
    MSG_SAVE_CDBG,

    // cDBG events
    MSG_ADD_DNODE,
    MSG_ADD_UNODE,
    MSG_DELETE_DNODE,
    MSG_DELETE_UNODE,
    MSG_ADD_EDGE,
    MSG_INCR_DNODE_COUNT
};


struct Event {
    Event(event_t msg_t)
        : msg_type(msg_t),
          msg(nullptr)
    {
    }

    Event(event_t msg_t, void * msg)
        : msg_type(msg_t),
          msg(msg)
    {
    }

    ~Event() {
        delete msg;
    }

    event_t msg_type;
    void * msg;
};



struct StreamingCompactorReport {
    uint64_t read_n;
    uint64_t n_full;
    uint64_t n_tips;
    uint64_t n_islands;
    uint64_t n_unknown;
    uint64_t n_trivial;
    uint64_t n_dnodes;
    uint64_t n_unodes;
    uint64_t n_updates;
    uint64_t n_tags;
    uint64_t n_unique;
    double   estimated_fp;
};


}
}

#endif
