/* cdbg/cdbg_types.hh -- cdbg node types
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_CDBG_TYPES_HH
#define BOINK_CDBG_TYPES_HH

#include <climits>
#include <memory>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include "boink/boink.hh"
#include "boink/hashing/alphabets.hh"
#include "boink/event_types.hh"

#define NULL_ID             ULLONG_MAX
#define UNITIG_START_ID     0

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

typedef uint64_t id_t;

enum node_meta_t {
    FULL,
    TIP,
    ISLAND,
    CIRCULAR,
    LOOP,
    TRIVIAL,
    DECISION
};


inline const char * node_meta_repr(node_meta_t meta) {
    switch(meta) {
        case FULL:
            return "FULL";
        case TIP:
            return "TIP";
        case ISLAND:
            return "ISLAND";
        case CIRCULAR:
            return "CIRCULAR";
        case TRIVIAL:
            return "TRIVIAL";
        case LOOP:
            return "LOOP";
        case DECISION:
            return "DECISION";
    }
}

enum cDBGFormat {
    GRAPHML,
    EDGELIST,
    FASTA,
    GFA1
};


inline const std::string cdbg_format_repr(cDBGFormat fmt) {
    switch(fmt) {
        case GRAPHML:
            return "graphml";
        case EDGELIST:
            return "edgelist";
        case FASTA:
            return "fasta";
        case GFA1:
            return "gfa1";
        default:
            return "FORMAT";
    }
}


enum update_meta_t {
    BUILD_UNODE,
    BUILD_DNODE,
    DELETE_UNODE,
    SPLIT_UNODE,
    EXTEND_UNODE,
    CLIP_UNODE,
    MERGE_UNODES
};



/* cDBG history events. These all encode the sequence of the new node,
 * the parents (if necessary), the children (if necessary), and related node
 * meta.
 */

struct HistoryNewEvent : public events::Event {
    HistoryNewEvent()
        : events::Event(events::MSG_HISTORY_NEW)
    {}

    std::string sequence;
    id_t id;
    node_meta_t meta;
};


struct HistoryMergeEvent : public events::Event {
    HistoryMergeEvent()
        : events::Event(events::MSG_HISTORY_MERGE)
    {}

    std::string sequence;
    id_t lparent, rparent, child;
    node_meta_t meta;
};


struct HistoryExtendEvent : public events::Event {
    HistoryExtendEvent()
        : events::Event(events::MSG_HISTORY_EXTEND)
    {}

    id_t id;
    std::string sequence;
    node_meta_t meta;
};


struct HistoryDeleteEvent : public events::Event {
    HistoryDeleteEvent()
        : events::Event(events::MSG_HISTORY_DELETE)
    {}

    id_t id;
    node_meta_t meta;
};


struct HistoryClipEvent : public events::Event {
    HistoryClipEvent()
        : events::Event(events::MSG_HISTORY_CLIP)
    {}

    id_t id;
    std::string sequence;
    node_meta_t meta;
};


struct HistorySplitEvent : public events::Event {
    HistorySplitEvent()
        : events::Event(events::MSG_HISTORY_SPLIT)
    {}

    id_t parent, lchild, rchild;
    node_meta_t lmeta, rmeta;
    std::string lsequence, rsequence;
};


struct HistorySplitCircularEvent : public events::Event {
    HistorySplitCircularEvent ()
        : events::Event(events::MSG_HISTORY_SPLIT_CIRCULAR)
    {}

    id_t id;
    std::string sequence;
    node_meta_t meta;
};



}
}

#undef pdebug
#endif
