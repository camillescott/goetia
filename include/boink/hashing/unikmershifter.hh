/**
 * (c) Camille Scott, 2019
 * File   : ukhshashshifter.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 10.01.2020
 */

#include "boink/hashing/ukhs.hh"
#include "boink/hashing/hashshifter.hh"

namespace boink::hashing {
    extern template class HashShifter<FwdUnikmerPolicy>;
    extern template class HashShifter<CanUnikmerPolicy>;
}