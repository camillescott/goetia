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
    
    template class signatures::UnikmerSignature<storage::BitStorage, hashing::HashModel<uint64_t>>;
    template class signatures::UnikmerSignature<storage::BitStorage, hashing::CanonicalModel<uint64_t>>;

    template class signatures::UnikmerSignature<storage::SparseppSetStorage, hashing::HashModel<uint64_t>>;
    template class signatures::UnikmerSignature<storage::SparseppSetStorage, hashing::CanonicalModel<uint64_t>>;

    template class signatures::UnikmerSignature<storage::ByteStorage, hashing::HashModel<uint64_t>>;
    template class signatures::UnikmerSignature<storage::ByteStorage, hashing::CanonicalModel<uint64_t>>;

    template class signatures::UnikmerSignature<storage::NibbleStorage, hashing::HashModel<uint64_t>>;
    template class signatures::UnikmerSignature<storage::NibbleStorage, hashing::CanonicalModel<uint64_t>>;

    template class signatures::UnikmerSignature<storage::QFStorage, hashing::HashModel<uint64_t>>;
    template class signatures::UnikmerSignature<storage::QFStorage, hashing::CanonicalModel<uint64_t>>;



}
