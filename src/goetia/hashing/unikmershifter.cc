/**
 * (c) Camille Scott, 2019
 * File   : unikmershifter.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 08.01.2020
 */

#include "goetia/hashing/ukhs.hh"
#include "goetia/hashing/unikmershifter.hh"
#include "goetia/hashing/hashshifter.hh"

namespace goetia {

    template class HashShifter<FwdUnikmerPolicy>;
    template HashShifter<FwdUnikmerPolicy>::HashShifter(uint16_t, uint16_t&&, std::shared_ptr<typename FwdUnikmerPolicy::ukhs_type>&&);
    template HashShifter<FwdUnikmerPolicy>::HashShifter(const std::string&, uint16_t, uint16_t&&, std::shared_ptr<typename FwdUnikmerPolicy::ukhs_type>&&);

    template class HashShifter<CanUnikmerPolicy>;   
    template HashShifter<CanUnikmerPolicy>::HashShifter(uint16_t, uint16_t&&, std::shared_ptr<typename CanUnikmerPolicy::ukhs_type>&&);
    template HashShifter<CanUnikmerPolicy>::HashShifter(const std::string&, uint16_t, uint16_t&&, std::shared_ptr<typename CanUnikmerPolicy::ukhs_type>&&);

}
