/**
 * (c) Camille Scott, 2019
 * File   : unikmershifter.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 08.01.2020
 */

#include "boink/hashing/ukhs.hh"
#include "boink/hashing/unikmershifter.hh"
#include "boink/hashing/hashshifter.hh"

namespace boink::hashing {
    template class HashShifter<FwdUnikmerPolicy>;
    template class HashShifter<CanUnikmerPolicy>;   
}