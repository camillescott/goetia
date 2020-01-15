/* kmerclient.hh: base class for things parametrized by K
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_SEQUENCES_EXCEPTIONS_HH
#define BOINK_SEQUENCES_EXCEPTIONS_HH

#include "boink/boink.hh"

namespace boink {

class InvalidCharacterException : public BoinkException {
public:
    explicit InvalidCharacterException(const std::string& msg = "Invalid character encountered.")
        : BoinkException(msg) { }
};


class SequenceLengthException: public BoinkException {
public:
    explicit SequenceLengthException(const std::string& msg = "Sequence was too short.")
        : BoinkException(msg) { }

};

}
#endif
