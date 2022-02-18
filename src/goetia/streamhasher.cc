/**
 * (c) Camille Scott, 2021
 * File   : streamhasher.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 17.02.2022
 */

#include "goetia/streamhasher.hh"

namespace goetia {
    template class StreamHasher<FwdLemireShifter>;
    template class StreamHasher<CanLemireShifter>;
}
