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

namespace goetia {

    template class dBGWalker<dBG<storage::BitStorage, hashing::FwdLemireShifter>>;
    template class dBGWalker<dBG<storage::BitStorage, hashing::CanLemireShifter>>;
    template class dBGWalker<dBG<storage::BitStorage, hashing::FwdUnikmerShifter>>;
    template class dBGWalker<dBG<storage::BitStorage, hashing::CanUnikmerShifter>>;

    template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::FwdLemireShifter>>;
    template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::CanLemireShifter>>;
    template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>>;
    template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>>;

    template class dBGWalker<dBG<storage::ByteStorage, hashing::FwdLemireShifter>>;
    template class dBGWalker<dBG<storage::ByteStorage, hashing::CanLemireShifter>>;
    template class dBGWalker<dBG<storage::ByteStorage, hashing::FwdUnikmerShifter>>;
    template class dBGWalker<dBG<storage::ByteStorage, hashing::CanUnikmerShifter>>;

    template class dBGWalker<dBG<storage::NibbleStorage, hashing::FwdLemireShifter>>;
    template class dBGWalker<dBG<storage::NibbleStorage, hashing::CanLemireShifter>>;
    template class dBGWalker<dBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>>;
    template class dBGWalker<dBG<storage::NibbleStorage, hashing::CanUnikmerShifter>>;

    template class dBGWalker<dBG<storage::QFStorage, hashing::FwdLemireShifter>>;
    template class dBGWalker<dBG<storage::QFStorage, hashing::CanLemireShifter>>;
    template class dBGWalker<dBG<storage::QFStorage, hashing::FwdUnikmerShifter>>;
    template class dBGWalker<dBG<storage::QFStorage, hashing::CanUnikmerShifter>>;

};
