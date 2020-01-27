/**
 * (c) Camille Scott, 2019
 * File   : hashshifter.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 21.01.2020
 */

#include <string>

#include "boink/hashing/hashshifter.hh"
#include "boink/hashing/rollinghashshifter.hh"


template class boink::hashing::HashShifter<boink::hashing::FwdLemirePolicy>;
template boink::hashing::HashShifter<boink::hashing::FwdLemirePolicy>::HashShifter(uint16_t);
template boink::hashing::HashShifter<boink::hashing::FwdLemirePolicy>::HashShifter(const std::string&, uint16_t);

template class boink::hashing::HashShifter<boink::hashing::CanLemirePolicy>;
template boink::hashing::HashShifter<boink::hashing::CanLemirePolicy>::HashShifter(uint16_t);
template boink::hashing::HashShifter<boink::hashing::CanLemirePolicy>::HashShifter(const std::string&, uint16_t);
