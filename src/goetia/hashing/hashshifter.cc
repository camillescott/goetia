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


template class goetia::HashShifter<goetia::FwdLemirePolicy>;
template goetia::HashShifter<goetia::FwdLemirePolicy>::HashShifter(uint16_t);
template goetia::HashShifter<goetia::FwdLemirePolicy>::HashShifter(const std::string&, uint16_t);

template class goetia::HashShifter<goetia::CanLemirePolicy>;
template goetia::HashShifter<goetia::CanLemirePolicy>::HashShifter(uint16_t);
template goetia::HashShifter<goetia::CanLemirePolicy>::HashShifter(const std::string&, uint16_t);
