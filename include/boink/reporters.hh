/* reporters.hh -- async reporters
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef REPORTERS_HH
#define REPORTERS_HH

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "boink/events.hh"
#include "boink/event_types.hh"

namespace boink {
namespace reporters {

using std::shared_ptr;

using boink::events::EventListener;
using boink::event_types::Event;
using boink::event_types::StreamingCompactorReport;

class Reporter : public EventListener {

protected:

    std::string _output_filename;
    std::ofstream _output_stream;

public:

    Reporter(const std::string& output_filename,
             const std::string thread_name) 
        : EventListener(thread_name),
          _output_filename(output_filename),
          _output_stream(_output_filename.c_str())
    {
    }

    Reporter(const std::string& output_filename)
        : Reporter(output_filename, "BaseReporter")
    {
    }

    virtual ~Reporter() {
        _output_stream.close();
    }
};


class StreamingCompactorReporter: public Reporter {

public:

    StreamingCompactorReporter(const std::string& output_filename)
        : Reporter(output_filename, "StreamingCompactorReporter") {
        
        this->msg_type_whitelist.insert(boink::event_types::MSG_WRITE_CDBG_STATS);

        _output_stream << "read_n,n_full,n_tips,n_islands,n_unknown"
                          ",n_trivial,n_dnodes,n_unodes,n_tags,"
                          "n_updates,n_unique,estimated_fp" << std::endl;
    }

    virtual void handle_msg(shared_ptr<Event> event) {
        if (event->msg_type == boink::event_types::MSG_WRITE_CDBG_STATS) {
            StreamingCompactorReport * report = static_cast<StreamingCompactorReport*>(event->msg);
            _output_stream << report->read_n << ","
                           << report->n_full << ","
                           << report->n_tips << ","
                           << report->n_islands << ","
                           << report->n_unknown << ","
                           << report->n_trivial << ","
                           << report->n_dnodes << ","
                           << report->n_unodes << ","
                           << report->n_tags << ","
                           << report->n_updates << ","
                           << report->n_unique << ","
                           << report->estimated_fp << ","
                           << std::endl;
        }
    }



};


}
}

#endif
