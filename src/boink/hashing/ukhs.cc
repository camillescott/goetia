/**
 * (c) Camille Scott, 2019
 * File   : ukhs.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.01.2020
 */

#include <string>

#include "boink/hashing/ukhs.hh"
#include "boink/hashing/rollinghashshifter.hh"

namespace boink::hashing {

    template class UKHS<FwdRollingShifter>;
    template class UKHS<CanRollingShifter>;
    template class UnikmerShifter<FwdRollingShifter>;
    template class UnikmerShifter<CanRollingShifter>;

    template<> template<>
    typename UnikmerShifter<FwdRollingShifter>::hash_type
    UnikmerShifter<FwdRollingShifter>::_hash_base<std::string::iterator>(std::string::iterator begin,
                                                                         std::string::iterator end);
    template<> template<>
    typename UnikmerShifter<CanRollingShifter>::hash_type
    UnikmerShifter<CanRollingShifter>::_hash_base<std::string::iterator>(std::string::iterator begin,
                                                                         std::string::iterator end);
}
