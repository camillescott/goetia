/**
 * (c) Camille Scott, 2019
 * File   : hashshifter.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 21.01.2020
 */

#include <string>

#include "goetia/hashing/hashshifter.hh"
#include "goetia/hashing/rollinghashshifter.hh"


template class goetia::hashing::HashShifter<goetia::hashing::FwdLemirePolicy>;
template goetia::hashing::HashShifter<goetia::hashing::FwdLemirePolicy>::HashShifter(uint16_t);
template goetia::hashing::HashShifter<goetia::hashing::FwdLemirePolicy>::HashShifter(const std::string&, uint16_t);

template class goetia::hashing::HashShifter<goetia::hashing::CanLemirePolicy>;
template goetia::hashing::HashShifter<goetia::hashing::CanLemirePolicy>::HashShifter(uint16_t);
template goetia::hashing::HashShifter<goetia::hashing::CanLemirePolicy>::HashShifter(const std::string&, uint16_t);
