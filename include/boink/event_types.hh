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

using std::string;

namespace boink {
namespace event_types {

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
    MSG_ADD_EDGE,
    MSG_INCR_DNODE_COUNT,
    MSG_EXTEND_UNODE,
    MSG_SPLIT_UNODE
};


struct Event {
    Event(event_t msg_t)
        : msg_type(msg_t)
    {
    }

    const event_t msg_type;
};


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


struct BuildDNodeEvent : public Event {
    BuildDNodeEvent()
        : Event(MSG_ADD_DNODE)
    {}
    hash_t hash;
    string kmer;
};


struct BuildUNodeEvent : public Event {
    BuildUNodeEvent()
        : Event(MSG_ADD_UNODE)
    {}
    HashVector tags;
    string sequence;
    junction_t left;
    junction_t right;
};


struct DeleteUNodeEvent : public Event {
    DeleteUNodeEvent()
        : Event(MSG_DELETE_UNODE)
    {}
    id_t node_id;
};


struct IncrDNodeEvent : public Event {
    IncrDNodeEvent()
        : Event(MSG_INCR_DNODE_COUNT)
    {}
    hash_t dnode;
};


struct ExtendUnodeEvent : public Event {
    ExtendUnodeEvent()
        : Event(MSG_EXTEND_UNODE)
    {}
    id_t node_id;
    direction_t from_dir;
    junction_t new_junction;
    string extra_sequence;
};


struct SplitUnodeEvent : public Event {
    SplitUnodeEvent()
        : Event(MSG_SPLIT_UNODE)
    {}
    id_t node_id;
    junction_t new_left;
    junction_t new_right;
    uint64_t split_position;
    hash_t dnode;
};

}
}

#endif
