/**
 * (c) Camille Scott, 2019
 * File   : hashextender.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 09.01.2020
 */

#include "boink/hashing/hashextender.hh"
#include "boink/hashing/shifter_types.hh"

namespace boink::hashing {

    template class DefaultExtensionPolicy<FwdRollingShifter>;
    template class DefaultExtensionPolicy<CanRollingShifter>;

    template class HashExtender<DefaultExtensionPolicy<FwdRollingShifter>>;
    template HashExtender<DefaultExtensionPolicy<FwdRollingShifter>>::HashExtender(uint16_t);
    template HashExtender<DefaultExtensionPolicy<FwdRollingShifter>>::HashExtender(const std::string&, uint16_t);

    template class HashExtender<DefaultExtensionPolicy<CanRollingShifter>>;
    template HashExtender<DefaultExtensionPolicy<CanRollingShifter>>::HashExtender(uint16_t);
    template HashExtender<DefaultExtensionPolicy<CanRollingShifter>>::HashExtender(const std::string&, uint16_t);

    template class HashExtender<FwdUnikmerShifter>;
    template class HashExtender<CanUnikmerShifter>;

    template class KmerIterator<FwdRollingExtender>;
    template class KmerIterator<CanRollingExtender>;

    template class KmerIterator<FwdUnikmerExtender>;
    template class KmerIterator<CanUnikmerExtender>;

}
