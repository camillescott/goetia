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

#include "boink/boink.hh"
#include "boink/cdbg.hh"
#include "boink/compactor.hh"
#include "boink/events.hh"
#include "boink/event_types.hh"

namespace boink {
namespace reporters {

using std::shared_ptr;

using boink::events::EventListener;
using boink::event_types::Event;
using boink::event_types::TimeIntervalEvent;

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


template <class GraphType>
class StreamingCompactorReporter: public Reporter {

protected:

    StreamingCompactor<GraphType> * compactor;

public:

    StreamingCompactorReporter(StreamingCompactor<GraphType> * compactor,
                               const std::string& output_filename)
        : Reporter(output_filename, "StreamingCompactorReporter"),
          compactor(compactor)
    {    
        this->msg_type_whitelist.insert(boink::event_types::MSG_TIME_INTERVAL);

        _output_stream << "read_n,n_full,n_tips,n_islands,n_unknown"
                          ",n_trivial,n_dnodes,n_unodes,n_tags,"
                          "n_updates,n_unique,estimated_fp" << std::endl;
    }

    virtual void handle_msg(shared_ptr<Event> event) {
        if (event->msg_type == boink::event_types::MSG_TIME_INTERVAL) {
            auto _event = static_cast<TimeIntervalEvent*>(event.get());
            if (_event->level == TimeIntervalEvent::FINE) {
                auto report = compactor->get_report();
                _output_stream << _event->t << ","
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
                               << report->estimated_fp 
                               << std::endl;
            }
        }
    }
};


class cDBGWriter : public Reporter {
protected:

    cDBG * cdbg;
    cDBGFormat format;

public:

    cDBGWriter(cDBG * cdbg,
               cDBGFormat format,
               const string& output_filename)
        : Reporter(output_filename, "cDBGWriter"),
          cdbg(cdbg),
          format(format)
    {
        this->msg_type_whitelist.insert(boink::event_types::MSG_TIME_INTERVAL);
    }

    virtual void handle_msg(shared_ptr<Event> event) {
        if (event->msg_type == boink::event_types::MSG_TIME_INTERVAL) {
            auto _event = static_cast<TimeIntervalEvent*>(event.get());
            if (_event->level == TimeIntervalEvent::MEDIUM) {
                cdbg->write(this->_output_stream, format);
            }
        }
    }

};

}
}

#endif
