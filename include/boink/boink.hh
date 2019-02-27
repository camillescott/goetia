/* boink.hh
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_MISC_HH
#define BOINK_MISC_HH

#include <exception>
#include <string>
#include <iostream>
#include <type_traits>
#include <iterator>

namespace boink {

#define _cerr(x) do { std::ostringstream stream; \
                      stream << x << std::endl; \
                      std::cerr << stream.str(); \
                    } while(0)
#define _cout(x) do { std::ostringstream stream; \
                      stream << x << std::endl; \
                      std::cout << stream.str(); \
                    } while(0)

enum direction_t {
    DIR_LEFT,
    DIR_RIGHT
};



template<typename _Ty1, typename _Ty2>
std::ostream& operator<<(std::ostream& _os, const std::pair<_Ty1, _Ty2>& _p) {
    _os << "(" << _p.first << ", " << _p.second << ")";
    return _os;
}

template <template<typename, typename> class Container, class T, class A>
std::ostream& operator<<(std::ostream& os, const Container<T,A>& v)
{
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i != v.size() - 1)
            os << ", ";
    }
    os << "]";
    return os;
}

class BoinkException : public std::exception {
public:
    explicit BoinkException(const std::string& msg = "Generic boink exception.")
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

} // boink


#endif
