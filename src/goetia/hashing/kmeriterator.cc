/**
 * (c) Camille Scott, 2019
 * File   : kmeriterator.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 29.08.2019
 */

#include "goetia/hashing/kmeriterator.hh"
#include "goetia/hashing/shifter_types.hh"

namespace goetia {

    template class KmerIterator<FwdLemireShifter>;
    template class KmerIterator<CanLemireShifter>;

    template class KmerIterator<FwdUnikmerShifter>;
    template class KmerIterator<CanUnikmerShifter>;

}
