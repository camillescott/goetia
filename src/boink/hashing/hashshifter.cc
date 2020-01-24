/**
 * (c) Camille Scott, 2019
 * File   : hashshifter.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 21.01.2020
 */


#include "boink/hashing/hashshifter.hh"
#include "boink/hashing/rollinghashshifter.hh"


template class boink::hashing::HashShifter<boink::hashing::FwdLemirePolicy>;
template class boink::hashing::HashShifter<boink::hashing::CanLemirePolicy>;