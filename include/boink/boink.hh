/* boink.hh
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_HH
#define BOINK_HH

#include <exception>
#include <string>
#include <iostream>
#include <vector>
#include <type_traits>
#include <iterator>
#include <set>

namespace boink {

class BoinkException : public std::exception {
public:
    explicit BoinkException(const std::string& msg = "Generic boink exception")
        : _msg(msg) { }

    virtual ~BoinkException() throw() { }
    virtual const char* what() const throw ()
    {
        return _msg.c_str();
    }

protected:
    const std::string _msg;
};

} // namespace boink


#endif
