/**
 * (c) Camille Scott, 2019
 * File   : goetia.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 13.01.2020
 */


#ifndef BOINK_MISC_HH
#define BOINK_MISC_HH

#include <cstdint>
#include <exception>
#include <string>
#include <iostream>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace goetia {

#define _cerr(x) do { std::ostringstream stream; \
                      stream << x << std::endl; \
                      std::cerr << stream.str(); \
                    } while(0)
#define _cout(x) do { std::ostringstream stream; \
                      stream << x << std::endl; \
                      std::cout << stream.str(); \
                    } while(0)

template <typename T>
std::string repr(const T& item) {
    std::ostringstream _os;
    _os << item;
    return _os.str();
}

template<typename _Ty1, typename _Ty2>
std::string repr(const std::pair<_Ty1, _Ty2>& _p) {
    std::ostringstream _os;
    _os << "(" << _p.first << ", " << _p.second << ")";
    return _os.str();
}

template <template<typename, typename> class Container, class T, class A>
std::string repr(const Container<T,A>& v) {
    std::ostringstream os;
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        os << repr(v[i]);
        if (i != v.size() - 1)
            os << ", ";
    }
    os << "]";
    return os.str();
}

template<template<typename, typename> class Container, class T, class A>
bool contains(const Container<T,A>& collection,
              T&                    item) {
    return collection.count(item);
}


template<typename T>
bool contains(std::vector<T> collection,
              T item) {
    for (auto& i : collection) {
        if (i == item) {
            return true;
        }
    }
    return false;
}


class BoinkException : public std::exception {
public:
    explicit BoinkException(const std::string& msg = "Generic goetia exception.")
        : _msg(msg) { }

    virtual ~BoinkException() throw() { }
    virtual const char* what() const throw ()
    {
        return _msg.c_str();
    }

protected:
    const std::string _msg;
};

class BoinkFileException: public BoinkException {
public:
    explicit BoinkFileException(const std::string& msg = "Error reading or writing file.")
        : BoinkException(msg) { }
};


class InvalidStream : public BoinkFileException 
{
public:
    InvalidStream()
        : BoinkFileException("Generic InvalidStream error") {}
    explicit InvalidStream(const std::string& msg)
        : BoinkFileException(msg) {}
};

class StreamReadError : public BoinkFileException
{
public:
    StreamReadError()
        : BoinkFileException("Generic StreamReadError error") {}
    explicit StreamReadError(const std::string& msg)
        : BoinkFileException(msg) {}
};


template <typename T>
struct fail : std::false_type 
{
};

} // goetia


#endif
