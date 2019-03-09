/* ukhs_signature_reporter.hh -- async reporters
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_UKHS_SIG_REPORTER_HH
#define BOINK_UKHS_SIG_REPORTER_HH

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "boink/boink.hh"
#include "boink/event_types.hh"
#include "boink/reporting/reporters.hh"
#include "boink/reporting/report_types.hh"
#include "boink/ukhs_signature.hh"

namespace boink {
namespace reporting {

class UKHSSignatureReporter : public SingleFileReporter {

private:

    std::shared_ptr<signatures::UKHSCountSignature> signature;

public:

    UKHSSignatureReporter(std::shared_ptr<signatures::UKHSCountSignature> signature,
                          const std::string&                              filename)
        : SingleFileReporter(filename, "UKHSSignatureReporter"),
          signature(signature)
    {
        _cerr(this->THREAD_NAME << " reporting at MEDIUM interval.");
        this->msg_type_whitelist.insert(events::MSG_TIME_INTERVAL);
    }

    virtual void handle_msg(std::shared_ptr<events::Event> event) {
         if (event->msg_type == events::MSG_TIME_INTERVAL) {
            auto _event = static_cast<events::TimeIntervalEvent*>(event.get());
            if (_event->level == events::TimeIntervalEvent::MEDIUM ||
                _event->level == events::TimeIntervalEvent::END) {
                
                _output_stream << _event->t;
                auto counts = signature->get_signature();
                for (auto& count : counts) {
                    _output_stream << ", " << count;
                }
                _output_stream << std::endl;
            }
        }       
    }
};

}
}

#endif
