/**
 * (c) Camille Scott, 2019
 * File   : exceptions.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 05.08.2019
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
