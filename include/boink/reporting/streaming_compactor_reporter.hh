/* streaming_compactor_reporter.hh -- async reporters
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_COMPACTOR_REPORTER_HH
#define BOINK_COMPACTOR_REPORTER_HH

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "boink/boink.hh"
#include "boink/event_types.hh"
#include "boink/cdbg/compactor.hh"
#include "boink/reporting/reporters.hh"
#include "boink/reporting/report_types.hh"

using namespace boink::cdbg;
using namespace boink::reporting::report_types;

namespace boink {
namespace reporting {

template <class GraphType>
class StreamingCompactorReporter: public SingleFileReporter {

protected:

    shared_ptr<StreamingCompactor<GraphType>> compactor;

public:

    StreamingCompactorReporter(shared_ptr<StreamingCompactor<GraphType>> compactor,
                               const std::string& output_filename)
        : SingleFileReporter(output_filename, "StreamingCompactorReporter"),
          compactor(compactor)
    {    
        _cerr(this->THREAD_NAME << " reporting at FINE interval.");

        this->msg_type_whitelist.insert(boink::event_types::MSG_TIME_INTERVAL);

        _output_stream << "read_n,n_full,n_tips,n_islands,n_trivial"
                          ",n_circular,n_loops,n_dnodes,n_unodes,n_tags,"
                          "n_updates,n_splits,n_merges,n_extends,n_clips"
                          "n_deletes,n_circular_merges,n_unique,estimated_fp" << std::endl;
    }

    virtual void handle_msg(shared_ptr<Event> event) {
        if (event->msg_type == boink::event_types::MSG_TIME_INTERVAL) {
            auto _event = static_cast<TimeIntervalEvent*>(event.get());
            if (_event->level == TimeIntervalEvent::FINE ||
                _event->level == TimeIntervalEvent::END) {
                auto report = compactor->get_report();
                _output_stream << _event->t << ","
                               << report.n_full << ","
                               << report.n_tips << ","
                               << report.n_islands << ","
                               << report.n_trivial << ","
                               << report.n_circular << ","
                               << report.n_loops << ","
                               << report.n_dnodes << ","
                               << report.n_unodes << ","
                               << report.n_tags << ","
                               << report.n_updates << ","
                               << report.n_splits << ","
                               << report.n_merges << ","
                               << report.n_extends << ","
                               << report.n_clips << ","
                               << report.n_deletes << ","
                               << report.n_circular_merges << ","
                               << report.n_unique << ","
                               << report.estimated_fp 
                               << std::endl;
            }
        }
    }
};


}
}

#endif
