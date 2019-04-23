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

#include <string>

#include "boink/boink.hh"


namespace boink {
namespace events {


/* All event types. 
 */

enum event_t {
    // EventListener control events
    MSG_EXIT_THREAD,
    MSG_TIMER,

    // Data events
    MSG_TIME_INTERVAL,

    // cDBG events
    MSG_ADD_DNODE,
    MSG_ADD_UNODE,
    MSG_DELETE_DNODE,
    MSG_DELETE_UNODE,
    MSG_SPLIT_UNODE,
    MSG_EXTEND_UNODE,
    MSG_MERGE_UNODES,
    MSG_INCR_DNODE_COUNT,
    MSG_HISTORY_NEW,
    MSG_HISTORY_SPLIT,
    MSG_HISTORY_MERGE,
    MSG_HISTORY_EXTEND,
    MSG_HISTORY_CLIP,
    MSG_HISTORY_SPLIT_CIRCULAR,
    MSG_HISTORY_DELETE
};

/* The base Event. All Events should inherit from
 * this, as well as define an event_t.
 *
 */

struct Event {
    Event(event_t msg_t)
        : msg_type(msg_t)
    {
    }

    const event_t msg_type;
};


/* Used for output timekeeping.
 *
 */

struct TimeIntervalEvent : public Event {
    TimeIntervalEvent()
        : Event(MSG_TIME_INTERVAL)
    {}

    enum interval_level_t {
        FINE,
        MEDIUM,
        COARSE,
        END
    };

    interval_level_t level;
    uint64_t t;
};

}
}

#endif
