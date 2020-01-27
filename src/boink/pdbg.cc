/**
 * (c) Camille Scott, 2019
 * File   : pdbg.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */

#include "boink/pdbg.hh"
#include "boink/storage/storage_types.hh"
#include "boink/hashing/shifter_types.hh"

namespace boink {

    template class PdBG<storage::BitStorage, hashing::FwdUnikmerShifter>;
    template class PdBG<storage::BitStorage, hashing::CanUnikmerShifter>;

    // template class PdBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>;
    // template class PdBG<storage::SparseppSetStorage, hashing::CanRollingShifter>;

    // template class PdBG<storage::ByteStorage, hashing::FwdRollingShifter>;
    // template class PdBG<storage::ByteStorage, hashing::CanRollingShifter>;

    // template class PdBG<storage::NibbleStorage, hashing::FwdRollingShifter>;
    // template class PdBG<storage::NibbleStorage, hashing::CanRollingShifter>;

    // template class PdBG<storage::QFStorage, hashing::FwdRollingShifter>;
    // template class PdBG<storage::QFStorage, hashing::CanRollingShifter>;

}
