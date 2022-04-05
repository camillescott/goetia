/**
 * (c) Camille Scott, 2021
 * File   : errors.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 04.04.2022
 */

#ifndef GOETIA_ERRORS_HH
#define GOETIA_ERRORS_HH

#include <exception>
#include <stdexcept>
#include <string>

struct GoetiaBaseException : public std::exception {};

struct InvalidSequence : public GoetiaBaseException {
    size_t read_number;
    const std::string sequence;

    InvalidSequence(size_t read_number, const std::string& sequence)
        : read_number(read_number), sequence(sequence)
    {}
};

struct InvalidRecord : public GoetiaBaseException {
    size_t read_number;
    const std::string sequence;

    InvalidRecord(size_t read_number, const std::string& sequence)
        : read_number(read_number), sequence(sequence)
    {}
};


struct InvalidRecordPair : public GoetiaBaseException {
    size_t pair_number;
    const std::string name_left, name_right;

    InvalidRecordPair(size_t number, const std::string& name_left, const std::string& name_right)
        : pair_number(number), name_left(name_left), name_right(name_right)
    {}
};


struct EndOfStream : public GoetiaBaseException {};

struct InvalidStream : public GoetiaBaseException {
    const std::string message;

    InvalidStream(const std::string& message)
        : message(message)
    {}
};

struct InvalidPairedStream : public GoetiaBaseException {};

struct StreamReadError : public GoetiaBaseException {
    size_t read_number;
    const std::string filename;

    StreamReadError(size_t read_number, const std::string& filename)
        : read_number(read_number), filename(filename)
    {}
};


struct SequenceTooShort : public GoetiaBaseException {
    const std::string sequence;
    SequenceTooShort(const std::string& sequence)
        : sequence(sequence)
    {}
};


struct UninitializedShifter : public GoetiaBaseException {};

struct NotImplemented : public GoetiaBaseException {};

struct InvalidPartition : public GoetiaBaseException {
    uint64_t partition;
    InvalidPartition(uint64_t partition)
        : partition(partition)
    {}
};

struct DeserializationError : public GoetiaBaseException {
    const std::string message;
    DeserializationError(const std::string& message)
        : message(message)
    {}
};

struct SerializationError : public GoetiaBaseException {
    const std::string message;
    SerializationError(const std::string& message)
        : message(message)
    {}
};


#endif
