/**
 * (c) Camille Scott, 2019
 * File   : exceptions.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 05.08.2019
 */

#ifndef GOETIA_SEQUENCES_EXCEPTIONS_HH
#define GOETIA_SEQUENCES_EXCEPTIONS_HH

#include "goetia/goetia.hh"

namespace goetia {

class InvalidCharacterException : public GoetiaException {
public:
    explicit InvalidCharacterException(const std::string& msg = "Invalid character encountered.")
        : GoetiaException(msg) { }
};


class SequenceLengthException: public GoetiaException {
public:
    explicit SequenceLengthException(const std::string& msg = "Sequence was too short.")
        : GoetiaException(msg) { }

};


class InvalidSequenceException: public GoetiaException {
public:
    explicit InvalidSequenceException(const std::string& msg = "Got invalid sequence (too short or malformed).")
        : GoetiaException(msg) { }

};

}
#endif
