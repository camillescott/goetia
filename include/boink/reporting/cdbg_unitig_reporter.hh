/* cdbg_unitig_reporter.hh -- async reporters
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_CDBG_UNITIG_REPORTER_HH
#define BOINK_CDBG_UNITIG_REPORTER_HH

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "utils/stringutils.h"

#include "boink/boink.hh"
#include "boink/cdbg/cdbg.hh"
#include "boink/cdbg/compactor.hh"
#include "boink/event_types.hh"
#include "boink/reporting/reporters.hh"
#include "boink/reporting/report_types.hh"

#include "sparsepp/spp.h"


namespace boink {
namespace reporting {


template <class GraphType>
class cDBGUnitigReporter : public SingleFileReporter {
private:

    std::shared_ptr<cdbg::cDBG<GraphType>>        cdbg;
    std::vector<size_t>                           bins;

public:

    cDBGUnitigReporter(std::shared_ptr<cdbg::cDBG<GraphType>> cdbg,
                       const std::string&                     filename,
                       std::vector<size_t>                    bins)
        : SingleFileReporter       (filename, "cDBGUnitigReporter"),
          cdbg                     (cdbg),
          bins                     (bins)
    {
        _cerr(this->THREAD_NAME << " reporting at MEDIUM interval.");
        this->msg_type_whitelist.insert(events::MSG_TIME_INTERVAL);
        
        _output_stream << "read_n";

        for (size_t bin = 0; bin < bins.size() - 1; bin++) {
            _output_stream << ", " << bins[bin] << "-" << bins[bin+1];
        }

        _output_stream << std::endl;
    }

    virtual void handle_msg(std::shared_ptr<events::Event> event) {
         if (event->msg_type == events::MSG_TIME_INTERVAL) {
            auto _event = static_cast<events::TimeIntervalEvent*>(event.get());
            if (_event->level == events::TimeIntervalEvent::MEDIUM ||
                _event->level == events::TimeIntervalEvent::END) {
                
                auto bin_sums = this->compute_bins();
                auto row      = utils::StringUtils::join(bin_sums, ", ");
                _output_stream << _event->t << ","
                               << row << ","
                               << std::endl;
            }
        }       
    }

    std::vector<size_t> compute_bins() {
        auto time_start = std::chrono::system_clock::now();
        auto lock       = cdbg->lock_nodes();
        _cerr("Summing unitig length bins...");

        std::vector<size_t> bin_sums(bins.size(), 0);

        for (auto it = cdbg->unodes_begin(); it != cdbg->unodes_end(); ++it) {
            auto seq_len = it->second->sequence.length();
            for (size_t bin_num = 0; bin_num < bins.size() - 1; bin_num++) {
                if (seq_len >= bins[bin_num] && seq_len < bins[bin_num+1]) {
                    bin_sums[bin_num] += seq_len;
                    break;
                }
            }
        }

        auto time_elapsed = std::chrono::system_clock::now() - time_start;
        _cerr("Finished summing unitig length bins. Elapsed time: " <<
              std::chrono::duration<double>(time_elapsed).count());

        return bin_sums;
    }
};


}
}
#endif
