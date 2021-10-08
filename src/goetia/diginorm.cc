/**
 * (c) Camille Scott, 2021
 * File   : diginorm.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 08.10.2021
 */

#include "goetia/diginorm.hh"
#include "goetia/storage/storage.hh"

namespace goetia {

    template class DiginormFilter<dBG<storage::NibbleStorage, hashing::FwdLemireShifter>>;
    template class DiginormFilter<dBG<storage::NibbleStorage, hashing::CanLemireShifter>>;

    template class DiginormFilter<dBG<storage::ByteStorage, hashing::FwdLemireShifter>>;
    template class DiginormFilter<dBG<storage::ByteStorage, hashing::CanLemireShifter>>;
}
