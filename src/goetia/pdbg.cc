/**
 * (c) Camille Scott, 2019
 * File   : pdbg.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */

#include "goetia/pdbg.hh"
#include "goetia/storage/storage_types.hh"
#include "goetia/hashing/shifter_types.hh"

namespace goetia {

    template class PdBG<BitStorage, FwdUnikmerShifter>;
    template class PdBG<BitStorage, CanUnikmerShifter>;

    template class PdBG<SparseppSetStorage, FwdUnikmerShifter>;
    template class PdBG<SparseppSetStorage, CanUnikmerShifter>;

    template class PdBG<PHMapStorage, FwdUnikmerShifter>;
    template class PdBG<PHMapStorage, CanUnikmerShifter>;

    template class PdBG<BTreeStorage, FwdUnikmerShifter>;
    template class PdBG<BTreeStorage, CanUnikmerShifter>;

    template class PdBG<ByteStorage, FwdUnikmerShifter>;
    template class PdBG<ByteStorage, CanUnikmerShifter>;

    template class PdBG<NibbleStorage, FwdUnikmerShifter>;
    template class PdBG<NibbleStorage, CanUnikmerShifter>;

    template class PdBG<QFStorage, FwdUnikmerShifter>;
    template class PdBG<QFStorage, CanUnikmerShifter>;

}
