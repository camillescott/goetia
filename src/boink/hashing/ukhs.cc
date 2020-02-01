/**
 * (c) Camille Scott, 2019
 * File   : ukhs.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.01.2020
 */

#include <string>

#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/hashshifter.hh"
#include "boink/hashing/ukhs.hh"


namespace boink::hashing {

    template class UKHS<FwdRollingShifter>;
    template class UKHS<CanRollingShifter>;
    template class UnikmerShifterPolicy<FwdLemirePolicy>;
    template class UnikmerShifterPolicy<CanLemirePolicy>;
    template class UnikmerLemirePolicy<Hash<uint64_t>, DNA_SIMPLE>;
    template class UnikmerLemirePolicy<Canonical<uint64_t>, DNA_SIMPLE>;

    template<> template<>
    typename UnikmerShifterPolicy<FwdLemirePolicy>::hash_type
    UnikmerShifterPolicy<FwdLemirePolicy>::hash_base_impl<std::string::iterator>(std::string::iterator begin,
                                                                                 std::string::iterator end);
    template<> template<>
    typename UnikmerShifterPolicy<CanLemirePolicy>::hash_type
    UnikmerShifterPolicy<CanLemirePolicy>::hash_base_impl<std::string::iterator>(std::string::iterator begin,
                                                                                 std::string::iterator end);
                                                                    
}
