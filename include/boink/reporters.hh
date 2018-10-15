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

class SingleFileReporter : public EventListener {

protected:

    std::string _output_filename;
    std::ofstream _output_stream;

public:

    SingleFileReporter(const std::string& output_filename,
             const std::string& thread_name) 
        : EventListener(thread_name),
          _output_filename(output_filename),
          _output_stream(_output_filename.c_str())
    {
    }

    SingleFileReporter(const std::string& output_filename)
        : SingleFileReporter(output_filename, "SingleFileReporter")
    {
    }

    virtual ~SingleFileReporter() {
        _output_stream.close();
    }
};


class MultiFileReporter : public EventListener {

protected:

    std::string _file_prefix;
    std::vector<std::string> _filenames;
    std::vector<std::ofstream> _streams;

public:

    MultiFileReporter(const std::string& prefix,
                      const std::string& thread_name)
        : EventListener(thread_name),
          _file_prefix(prefix)
    {
    }

    MultiFileReporter(const std::string& prefix)
        : MultiFileReporter(prefix, "MultiFileReporter")
    {
    }

    std::ofstream& current_stream() {
        return _streams.back();
    }

    std::string& current_filename() {
        return _filenames.back();
    }

    std::ofstream& next_stream(uint64_t start_time, 
                               const std::string& suffix) {
        std::ostringstream name;
        name << _file_prefix << "." << start_time << "."
             << suffix;
        std::string filename = name.str();

        if (_streams.size() > 0) {
            current_stream().close();   
        }
        _streams.push_back(std::move(std::ofstream(filename.c_str())));
        _filenames.push_back(filename);
        return _streams.back();
    }

    virtual ~MultiFileReporter() {
        current_stream().close();
    }

};


template <class GraphType>
class StreamingCompactorReporter: public SingleFileReporter {

protected:

    StreamingCompactor<GraphType> * compactor;

public:

    StreamingCompactorReporter(StreamingCompactor<GraphType> * compactor,
                               const std::string& output_filename)
        : SingleFileReporter(output_filename, "StreamingCompactorReporter"),
          compactor(compactor)
    {    
        this->msg_type_whitelist.insert(boink::event_types::MSG_TIME_INTERVAL);

        _output_stream << "read_n,n_full,n_tips,n_islands,n_trivial"
                          ",n_circular,n_dnodes,n_unodes,n_tags,"
                          "n_updates,n_unique,estimated_fp" << std::endl;
    }

    virtual void handle_msg(shared_ptr<Event> event) {
        if (event->msg_type == boink::event_types::MSG_TIME_INTERVAL) {
            auto _event = static_cast<TimeIntervalEvent*>(event.get());
            if (_event->level == TimeIntervalEvent::FINE ||
                _event->level == TimeIntervalEvent::END) {
                auto report = compactor->get_report();
                _output_stream << _event->t << ","
                               << report->n_full << ","
                               << report->n_tips << ","
                               << report->n_islands << ","
                               << report->n_trivial << ","
                               << report->n_circular << ","
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


template <class GraphType>
class cDBGWriter : public MultiFileReporter {
protected:

    cDBG<GraphType> * cdbg;
    cDBGFormat format;

public:

    cDBGWriter(cDBG<GraphType> * cdbg,
               cDBGFormat format,
               const string& output_prefix)
        : MultiFileReporter(output_prefix,
                            "cDBGWriter[" + cdbg_format_repr(format) + "]"),
          cdbg(cdbg),
          format(format)
    {
        this->msg_type_whitelist.insert(boink::event_types::MSG_TIME_INTERVAL);
    }

    virtual void handle_msg(shared_ptr<Event> event) {
        if (event->msg_type == boink::event_types::MSG_TIME_INTERVAL) {
            auto _event = static_cast<TimeIntervalEvent*>(event.get());
            if (_event->level == TimeIntervalEvent::MEDIUM ||
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
