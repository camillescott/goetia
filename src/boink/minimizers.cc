/**
 * (c) Camille Scott, 2019
 * File   : minimizers.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.01.2020
 */

#include "boink/minimizers.hh"
#include "boink/hashing/rollinghashshifter.hh"

namespace boink {

template class InteriorMinimizer<uint64_t>;
template class WKMinimizer<hashing::FwdRollingShifter>;
template class WKMinimizer<hashing::CanRollingShifter>;

}
