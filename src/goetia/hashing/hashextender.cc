/**
 * (c) Camille Scott, 2019
 * File   : hashextender.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 09.01.2020
 */

#include "goetia/hashing/hashextender.hh"
#include "goetia/hashing/shifter_types.hh"

namespace goetia::hashing {

    template class DefaultExtensionPolicy<FwdLemireShifter>;
    template class DefaultExtensionPolicy<CanLemireShifter>;

    template class HashExtender<DefaultExtensionPolicy<FwdLemireShifter>>;
    template HashExtender<DefaultExtensionPolicy<FwdLemireShifter>>::HashExtender(uint16_t);
    template HashExtender<DefaultExtensionPolicy<FwdLemireShifter>>::HashExtender(const std::string&, uint16_t);

    template class HashExtender<DefaultExtensionPolicy<CanLemireShifter>>;
    template HashExtender<DefaultExtensionPolicy<CanLemireShifter>>::HashExtender(uint16_t);
    template HashExtender<DefaultExtensionPolicy<CanLemireShifter>>::HashExtender(const std::string&, uint16_t);

    template class HashExtender<FwdUnikmerShifter>;
    template class HashExtender<CanUnikmerShifter>;

    template class KmerIterator<FwdRollingExtender>;
    template class KmerIterator<CanRollingExtender>;

    template class KmerIterator<FwdUnikmerExtender>;
    template class KmerIterator<CanUnikmerExtender>;

}
