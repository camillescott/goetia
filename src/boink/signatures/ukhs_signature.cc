/**
 * (c) Camille Scott, 2019
 * File   : ukhs_signature.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 31.08.2019
 */


#include "boink/signatures/ukhs_signature.hh"

#include "boink/storage/storage_types.hh"
#include "boink/hashing/rollinghashshifter.hh"


namespace boink {

    template class signatures::UnikmerSignature<storage::BitStorage, hashing::FwdRollingShifter>;
    template class signatures::UnikmerSignature<storage::BitStorage, hashing::CanRollingShifter>;

    template class signatures::UnikmerSignature<storage::SparseppSetStorage, hashing::FwdRollingShifter>;
    template class signatures::UnikmerSignature<storage::SparseppSetStorage, hashing::CanRollingShifter>;

    template class signatures::UnikmerSignature<storage::ByteStorage, hashing::FwdRollingShifter>;
    template class signatures::UnikmerSignature<storage::ByteStorage, hashing::CanRollingShifter>;

    template class signatures::UnikmerSignature<storage::NibbleStorage, hashing::FwdRollingShifter>;
    template class signatures::UnikmerSignature<storage::NibbleStorage, hashing::CanRollingShifter>;

    template class signatures::UnikmerSignature<storage::QFStorage, hashing::FwdRollingShifter>;
    template class signatures::UnikmerSignature<storage::QFStorage, hashing::CanRollingShifter>;

}
