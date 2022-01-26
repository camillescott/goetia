/**
 * (c) Camille Scott, 2019
 * File   : assembly.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */

#include "goetia/dbg.hh"
#include "goetia/storage/storage_types.hh"
#include "goetia/hashing/rollinghashshifter.hh"
#include "goetia/hashing/ukhs.hh"
#include "goetia/traversal/unitig_walker.hh"

namespace goetia {

    template class UnitigWalker<dBG<BitStorage, FwdLemireShifter>>;
    template class UnitigWalker<dBG<BitStorage, CanLemireShifter>>;
    template class UnitigWalker<dBG<BitStorage, FwdUnikmerShifter>>;
    template class UnitigWalker<dBG<BitStorage, CanUnikmerShifter>>;

    template class UnitigWalker<dBG<SparseppSetStorage, FwdLemireShifter>>;
    template class UnitigWalker<dBG<SparseppSetStorage, CanLemireShifter>>;
    template class UnitigWalker<dBG<SparseppSetStorage, FwdUnikmerShifter>>;
    template class UnitigWalker<dBG<SparseppSetStorage, CanUnikmerShifter>>;

    template class UnitigWalker<dBG<ByteStorage, FwdLemireShifter>>;
    template class UnitigWalker<dBG<ByteStorage, CanLemireShifter>>;
    template class UnitigWalker<dBG<ByteStorage, FwdUnikmerShifter>>;
    template class UnitigWalker<dBG<ByteStorage, CanUnikmerShifter>>;

    template class UnitigWalker<dBG<NibbleStorage, FwdLemireShifter>>;
    template class UnitigWalker<dBG<NibbleStorage, CanLemireShifter>>;
    template class UnitigWalker<dBG<NibbleStorage, FwdUnikmerShifter>>;
    template class UnitigWalker<dBG<NibbleStorage, CanUnikmerShifter>>;

    template class UnitigWalker<dBG<QFStorage, FwdLemireShifter>>;
    template class UnitigWalker<dBG<QFStorage, CanLemireShifter>>;
    template class UnitigWalker<dBG<QFStorage, FwdUnikmerShifter>>;
    template class UnitigWalker<dBG<QFStorage, CanUnikmerShifter>>;

    template class goetia::UnitigWalker<goetia::dBG<goetia::PHMapStorage, goetia::FwdLemireShifter>>;
    template class goetia::UnitigWalker<goetia::dBG<goetia::PHMapStorage, goetia::CanLemireShifter>>;
    template class goetia::UnitigWalker<goetia::dBG<goetia::PHMapStorage, goetia::FwdUnikmerShifter>>;
    template class goetia::UnitigWalker<goetia::dBG<goetia::PHMapStorage, goetia::CanUnikmerShifter>>;

    template class goetia::UnitigWalker<goetia::dBG<goetia::BTreeStorage, goetia::FwdLemireShifter>>;
    template class goetia::UnitigWalker<goetia::dBG<goetia::BTreeStorage, goetia::CanLemireShifter>>;
    template class goetia::UnitigWalker<goetia::dBG<goetia::BTreeStorage, goetia::FwdUnikmerShifter>>;
    template class goetia::UnitigWalker<goetia::dBG<goetia::BTreeStorage, goetia::CanUnikmerShifter>>;



};
