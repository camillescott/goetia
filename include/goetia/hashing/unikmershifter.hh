/**
 * (c) Camille Scott, 2019
 * File   : ukhshashshifter.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 10.01.2020
 */

#include "goetia/hashing/ukhs.hh"
#include "goetia/hashing/hashshifter.hh"

namespace goetia {
    extern template class HashShifter<FwdUnikmerPolicy>;
    extern template class HashShifter<CanUnikmerPolicy>;

    typedef HashShifter<FwdUnikmerPolicy> FwdUnikmerShifter;
    typedef HashShifter<CanUnikmerPolicy> CanUnikmerShifter;
}
