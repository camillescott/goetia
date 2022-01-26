/**
 * (c) Camille Scott, 2019
 * File   : ukhs.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.01.2020
 */

#include <string>

#include "goetia/hashing/rollinghashshifter.hh"
#include "goetia/hashing/hashshifter.hh"
#include "goetia/hashing/ukhs.hh"


namespace goetia {

    template class UKHS<FwdLemireShifter>;
    template class UKHS<CanLemireShifter>;
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
