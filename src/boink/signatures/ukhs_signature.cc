/**
 * (c) Camille Scott, 2019
 * File   : ukhs_signature.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 31.08.2019
 */


#include "boink/signatures/ukhs_signature.hh"

#include "boink/storage/storage_types.hh"
#include "boink/hashing/canonical.hh"


namespace boink {
    
    template class signatures::UnikmerSignature<storage::BitStorage, hashing::Hash<uint64_t>>;
    template class signatures::UnikmerSignature<storage::BitStorage, hashing::Canonical<uint64_t>>;

    template class signatures::UnikmerSignature<storage::SparseppSetStorage, hashing::Hash<uint64_t>>;
    template class signatures::UnikmerSignature<storage::SparseppSetStorage, hashing::Canonical<uint64_t>>;

    template class signatures::UnikmerSignature<storage::ByteStorage, hashing::Hash<uint64_t>>;
    template class signatures::UnikmerSignature<storage::ByteStorage, hashing::Canonical<uint64_t>>;

    template class signatures::UnikmerSignature<storage::NibbleStorage, hashing::Hash<uint64_t>>;
    template class signatures::UnikmerSignature<storage::NibbleStorage, hashing::Canonical<uint64_t>>;

    template class signatures::UnikmerSignature<storage::QFStorage, hashing::Hash<uint64_t>>;
    template class signatures::UnikmerSignature<storage::QFStorage, hashing::Canonical<uint64_t>>;



}
