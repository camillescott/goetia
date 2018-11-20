/* cdbg_history_reporter.hh -- async reporters
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_CDBG_COMPONENT_REPORTER_HH
#define BOINK_CDBG_COMPONENT_REPORTER_HH

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "boink/boink.hh"
#include "boink/cdbg.hh"
#include "boink/compactor.hh"
#include "boink/event_types.hh"
#include "boink/reporting/reporters.hh"
#include "boink/reporting/report_types.hh"

#include <sparsepp/sparsepp/spp.h>

namespace boink {
namespace reporting {

class cDBGComponentReporter : public SingleFileReporter {
private:

    id_t _component_id_counter;

    class WeakComponent {
        public:

            id_t component_id;
            spp::sparse_hash_set<id_t> nodes;

            WeakComponent(id_t id)
                : component_id(id)
            {
            }

    };

    //spp::sparse_hash_map<id_t, unique_ptr<WeakComponent>> component_owner;
    //spp::sparse_hash_map<id_t, *WeakComponent> node_component_map;


public:

    cDBGComponentReporter(const std::string& filename)
        : SingleFileReporter(filename, "cDBGComponentReporter")
    {
        _cerr(this->THREAD_NAME << " reporting continuously.");

        this->msg_type_whitelist.insert(boink::event_types::MSG_HISTORY_NEW);
        this->msg_type_whitelist.insert(boink::event_types::MSG_HISTORY_SPLIT);
        this->msg_type_whitelist.insert(boink::event_types::MSG_HISTORY_SPLIT_CIRCULAR);
        this->msg_type_whitelist.insert(boink::event_types::MSG_HISTORY_MERGE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_HISTORY_EXTEND);
        this->msg_type_whitelist.insert(boink::event_types::MSG_HISTORY_CLIP);
        this->msg_type_whitelist.insert(boink::event_types::MSG_HISTORY_DELETE);

    }

    virtual void handle_msg(shared_ptr<Event> event) {
        if (event->msg_type == boink::event_types::MSG_HISTORY_NEW) {
            auto _event = static_cast<HistoryNewEvent*>(event.get());

        } else if (event->msg_type == boink::event_types::MSG_HISTORY_SPLIT) {
            auto _event = static_cast<HistorySplitEvent*>(event.get());

        } else if (event->msg_type == boink::event_types::MSG_HISTORY_MERGE) {
            auto _event = static_cast<HistoryMergeEvent*>(event.get());

        } else if (event->msg_type == boink::event_types::MSG_HISTORY_EXTEND) {
            auto _event = static_cast<HistoryExtendEvent*>(event.get());

        } else if (event->msg_type == boink::event_types::MSG_HISTORY_CLIP) {
            auto _event = static_cast<HistoryClipEvent*>(event.get());

        } else if (event->msg_type == boink::event_types::MSG_HISTORY_SPLIT_CIRCULAR) {
            auto _event = static_cast<HistorySplitCircularEvent*>(event.get());

        }
    }
};

}
}

#endif
