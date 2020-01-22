/**
 * (c) Camille Scott, 2019
 * File   : hashextender.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 09.01.2020
 */

#include "boink/hashing/hashextender.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/ukhs.hh"

namespace boink::hashing {

    template class HashExtender<FwdRollingShifter>;
    template class HashExtender<CanRollingShifter>;

    template class HashExtender<FwdUnikmerShifter>;
    template class HashExtender<CanUnikmerShifter>;

    template class KmerIterator<HashExtender<FwdRollingShifter>>;
    template class KmerIterator<HashExtender<CanRollingShifter>>;

    template class KmerIterator<HashExtender<FwdUnikmerShifter>>;
    template class KmerIterator<HashExtender<CanUnikmerShifter>>;

}
