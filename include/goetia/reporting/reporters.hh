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

#include "goetia/events.hh"

namespace goetia {

namespace reporting {

class SingleFileReporter : public events::EventListener {

protected:

    std::string _output_filename;
    std::ofstream _output_stream;

public:

    SingleFileReporter(const std::string& output_filename,
                       const std::string& thread_name);

    SingleFileReporter(const std::string& output_filename);

    virtual ~SingleFileReporter();
};


class MultiFileReporter : public events::EventListener {

protected:

    std::string _file_prefix;
    std::vector<std::string> _filenames;
    std::vector<std::ofstream> _streams;

public:

    MultiFileReporter(const std::string& prefix,
                      const std::string& thread_name);

    MultiFileReporter(const std::string& prefix);

    std::ofstream& current_stream();

    std::string& current_filename();

    std::ofstream& next_stream(uint64_t start_time, 
                               const std::string& suffix);

    virtual ~MultiFileReporter();

};

}
}

#endif
