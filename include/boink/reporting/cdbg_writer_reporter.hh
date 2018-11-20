/* cdbg_writer_reporter.hh -- async reporters
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_CDBG_WRITER_REPORTER_HH
#define BOINK_CDBG_WRITER_REPORTER_HH

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

template <class GraphType>
class cDBGWriter : public MultiFileReporter {
protected:

    shared_ptr<cDBG<GraphType>> cdbg;
    cDBGFormat format;

public:

    cDBGWriter(shared_ptr<cDBG<GraphType>> cdbg,
               cDBGFormat format,
               const string& output_prefix)
        : MultiFileReporter(output_prefix,
                            "cDBGWriter[" + cdbg_format_repr(format) + "]"),
          cdbg(cdbg),
          format(format)
    {
        _cerr(this->THREAD_NAME << " reporting at COARSE interval.");

        this->msg_type_whitelist.insert(boink::event_types::MSG_TIME_INTERVAL);
    }

    virtual void handle_msg(shared_ptr<Event> event) {
        if (event->msg_type == boink::event_types::MSG_TIME_INTERVAL) {
            auto _event = static_cast<TimeIntervalEvent*>(event.get());
            if (_event->level == TimeIntervalEvent::COARSE ||
                _event->level == TimeIntervalEvent::END) {

                std::ofstream& stream = this->next_stream(_event->t,
                                                          cdbg_format_repr(format));
                std::string&   filename = this->current_filename();

                _cerr(this->THREAD_NAME << ", t=" << _event->t <<
                      ": write cDBG to " << filename);
                cdbg->write(stream, format);
            }
        }
    }
};


}
}

#endif
