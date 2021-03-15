/**
 * (c) Camille Scott, 2019
 * File   : ukhs_signature.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 31.08.2019
 */


#include "goetia/sketches/unikmer_sketch.hh"

#include "goetia/storage/storage_types.hh"
#include "goetia/sketches/hllcounter.hh"
#include "goetia/hashing/canonical.hh"


namespace goetia {
    
    template class sketches::UnikmerSketch<storage::BitStorage, hashing::Hash<uint64_t>>;
    template class sketches::UnikmerSketch<storage::BitStorage, hashing::Canonical<uint64_t>>;

    template class sketches::UnikmerSketch<storage::SparseppSetStorage, hashing::Hash<uint64_t>>;
    template class sketches::UnikmerSketch<storage::SparseppSetStorage, hashing::Canonical<uint64_t>>;

    template class sketches::UnikmerSketch<storage::PHMapStorage, hashing::Hash<uint64_t>>;
    template class sketches::UnikmerSketch<storage::PHMapStorage, hashing::Canonical<uint64_t>>;

    template class sketches::UnikmerSketch<storage::BTreeStorage, hashing::Hash<uint64_t>>;
    template class sketches::UnikmerSketch<storage::BTreeStorage, hashing::Canonical<uint64_t>>;

    template class sketches::UnikmerSketch<storage::ByteStorage, hashing::Hash<uint64_t>>;
    template class sketches::UnikmerSketch<storage::ByteStorage, hashing::Canonical<uint64_t>>;

    template class sketches::UnikmerSketch<storage::NibbleStorage, hashing::Hash<uint64_t>>;
    template class sketches::UnikmerSketch<storage::NibbleStorage, hashing::Canonical<uint64_t>>;

    template class sketches::UnikmerSketch<storage::QFStorage, hashing::Hash<uint64_t>>;
    template class sketches::UnikmerSketch<storage::QFStorage, hashing::Canonical<uint64_t>>;

    template class sketches::UnikmerSketch<storage::HLLStorage, hashing::Hash<uint64_t>>;
    template class sketches::UnikmerSketch<storage::HLLStorage, hashing::Canonical<uint64_t>>;
}
