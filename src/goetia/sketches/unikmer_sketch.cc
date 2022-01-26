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
    
    template class UnikmerSketch<BitStorage, Hash<uint64_t>>;
    template class UnikmerSketch<BitStorage, Canonical<uint64_t>>;

    template class UnikmerSketch<SparseppSetStorage, Hash<uint64_t>>;
    template class UnikmerSketch<SparseppSetStorage, Canonical<uint64_t>>;

    template class UnikmerSketch<PHMapStorage, Hash<uint64_t>>;
    template class UnikmerSketch<PHMapStorage, Canonical<uint64_t>>;

    template class UnikmerSketch<BTreeStorage, Hash<uint64_t>>;
    template class UnikmerSketch<BTreeStorage, Canonical<uint64_t>>;

    template class UnikmerSketch<ByteStorage, Hash<uint64_t>>;
    template class UnikmerSketch<ByteStorage, Canonical<uint64_t>>;

    template class UnikmerSketch<NibbleStorage, Hash<uint64_t>>;
    template class UnikmerSketch<NibbleStorage, Canonical<uint64_t>>;

    template class UnikmerSketch<QFStorage, Hash<uint64_t>>;
    template class UnikmerSketch<QFStorage, Canonical<uint64_t>>;

    template class UnikmerSketch<HLLStorage, Hash<uint64_t>>;
    template class UnikmerSketch<HLLStorage, Canonical<uint64_t>>;
}
