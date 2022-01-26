/**
 * (c) Camille Scott, 2019
 * File   : minimizers.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.01.2020
 */

#include "goetia/minimizers.hh"
#include "goetia/hashing/rollinghashshifter.hh"

namespace goetia {

template class InteriorMinimizer<uint64_t>;
template class WKMinimizer<FwdLemireShifter>;
template class WKMinimizer<CanLemireShifter>;

}
