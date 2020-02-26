/* goetia.hh
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "goetia/goetia.hh"
#include "goetia/reporting/reporters.hh"

namespace goetia {
namespace reporting {

SingleFileReporter::SingleFileReporter(const std::string& output_filename,
                   const std::string& thread_name) 
    : EventListener(thread_name),
      _output_filename(output_filename),
      _output_stream(_output_filename.c_str())
{
}

SingleFileReporter::SingleFileReporter(const std::string& output_filename)
    : SingleFileReporter(output_filename, std::string("SingleFileReporter"))
{
}

SingleFileReporter::~SingleFileReporter() {
    _output_stream.close();

}

MultiFileReporter::MultiFileReporter(const std::string& prefix,
                                     const std::string& thread_name)
    : EventListener(thread_name),
      _file_prefix(prefix)
{
}

MultiFileReporter::MultiFileReporter(const std::string& prefix)
    : MultiFileReporter(prefix, std::string("MultiFileReporter"))
{
}

std::ofstream& MultiFileReporter::current_stream() {
    return _streams.back();
}

std::string& MultiFileReporter::current_filename() {
    return _filenames.back();
}

std::ofstream& MultiFileReporter::next_stream(uint64_t start_time, 
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

MultiFileReporter::~MultiFileReporter() {
    current_stream().close();
}

}
}
