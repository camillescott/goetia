/**
 * (c) Camille Scott, 2019
 * File   : ukhshashshifter.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 08.01.2020
 */

#include "boink/hashing/ukhshashshifter.hh"

namespace boink {

namespace hashing {

template class UnikmerShifter<FwdRollingShifter>;
template class UnikmerShifter<CanRollingShifter>;

}

template<> const std::string Tagged<hashing::FwdUnikmerShifter>::NAME = "Boink::UnikmerShifter::CyclicHash<uint64_t>::NonRandom_Fwd";
template<> const std::string Tagged<hashing::CanUnikmerShifter>::NAME = "Boink::UnikmerShifter::CyclicHash<uint64_t>::NonRandom_Can";


}
