/**
 * (c) Camille Scott, 2019
 * File   : assembly.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */

#include "boink/dbg.hh"
#include "boink/storage/storage_types.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/ukhs.hh"

namespace boink {

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
