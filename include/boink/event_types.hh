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
    MSG_SPLIT_UNODE,
    MSG_EXTEND_UNODE,
    MSG_MERGE_UNODES,
    MSG_INCR_DNODE_COUNT
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
    bool has_left;
    hash_t left_end;
    bool has_right;
    hash_t right_end;
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


struct ExtendUNodeEvent : public Event {
    ExtendUNodeEvent()
        : Event(MSG_EXTEND_UNODE)
    {}
    hash_t tip_hash;
    string new_sequence;
    direction_t from;
};


struct MergeUNodeEvent : public Event {
    MergeUNodeEvent()
        : Event(MSG_MERGE_UNODES)
    {}
    hash_t left_tip;
    hash_t right_tip;
    string new_sequence;
};


struct SplitUNodeEvent : public Event {
    SplitUNodeEvent()
        : Event(MSG_SPLIT_UNODE)
    {}
    id_t unode_id;
    size_t split_at;
};

}
}

#endif
