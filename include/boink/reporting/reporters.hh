/* reporters.hh -- async reporter base classes.
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_REPORTERS_HH
#define BOINK_REPORTERS_HH

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "boink/boink.hh"
#include "boink/events.hh"
#include "boink/event_types.hh"

namespace boink {
namespace reporting {

using std::shared_ptr;

using boink::events::EventListener;
using boink::events::Event;
using boink::events::TimeIntervalEvent;

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


}
}

#endif
