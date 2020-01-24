/**
 * (c) Camille Scott, 2019
 * File   : kmeriterator.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 29.08.2019
 */

#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/shifter_types.hh"

namespace boink::hashing {

    template class KmerIterator<FwdRollingShifter>;
    template class KmerIterator<CanRollingShifter>;

    template class KmerIterator<FwdUnikmerShifter>;
    template class KmerIterator<CanUnikmerShifter>;

}
