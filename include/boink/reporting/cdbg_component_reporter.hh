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
#include "boink/cdbg/cdbg.hh"
#include "boink/cdbg/compactor.hh"
#include "boink/event_types.hh"
#include "boink/reporting/reporters.hh"
#include "boink/reporting/report_types.hh"

#include <sparsepp/sparsepp/spp.h>

using namespace boink::cdbg;

namespace boink {
namespace reporting {


template <class GraphType>
class cDBGComponentReporter : public SingleFileReporter {
private:

    std::shared_ptr<cDBG<GraphType>> cdbg;
    id_t _component_id_counter;
    spp::sparse_hash_map<id_t, uint64_t> component_size_map;
    uint64_t max_component;
    uint64_t min_component;

public:

    cDBGComponentReporter(std::shared_ptr<cDBG<GraphType>> cdbg,
                          const std::string& filename)
        : SingleFileReporter(filename, "cDBGComponentReporter"),
          cdbg(cdbg),
          min_component(ULLONG_MAX),
          max_component(0)
    {
        _cerr(this->THREAD_NAME << " reporting at MEDIUM interval.");
        this->msg_type_whitelist.insert(boink::event_types::MSG_TIME_INTERVAL);

        _output_stream << "read_n,n_components,max_component,min_component" << std::endl;
    }

    virtual void handle_msg(std::shared_ptr<Event> event) {
         if (event->msg_type == boink::event_types::MSG_TIME_INTERVAL) {
            auto _event = static_cast<TimeIntervalEvent*>(event.get());
            if (_event->level == TimeIntervalEvent::MEDIUM ||
                _event->level == TimeIntervalEvent::END) {
                
                this->recompute_components();
                _output_stream << _event->t << ","
                               << component_size_map.size() << ","
                               << max_component << ","
                               << min_component << ","
                               << std::endl;
            }
        }       
    }

    void recompute_components() {
        auto lock = cdbg->lock_nodes();

        spp::sparse_hash_set<id_t> seen;
        component_size_map.clear();

        for (auto unitig_it = cdbg->unodes_begin(); unitig_it != cdbg->unodes_end(); unitig_it++) {
            auto root = unitig_it->second.get();
            if (!seen.count(root->node_id)) {
                auto component_nodes = cdbg->traverse_breadth_first(root);
                id_t component_id = (root->component_id == NULL_ID) ? ++_component_id_counter : root->component_id;
                for (auto node : component_nodes) {
                    node->component_id = component_id;
                    seen.insert(node->node_id);
                }
                component_size_map[component_id] = component_nodes.size();
                max_component = (component_nodes.size() > max_component) ? component_nodes.size() : max_component;
                min_component = (component_nodes.size() < min_component) ? component_nodes.size() : min_component;
            }
        }
    }
};

}
}

#endif
