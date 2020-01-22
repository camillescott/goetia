/**
 * (c) Camille Scott, 2019
 * File   : ukhs.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.01.2020
 */


#include "boink/hashing/ukhs.hh"
#include "boink/hashing/rollinghashshifter.hh"

namespace boink::hashing {

    template class UKHS<FwdRollingShifter>;
    template class UKHS<CanRollingShifter>;
    template class UnikmerShifter<FwdRollingShifter>;
    template class UnikmerShifter<CanRollingShifter>;

}
