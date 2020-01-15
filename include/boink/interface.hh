/**
 * (c) Camille Scott, 2019
 * File   : interface.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 25.09.2019
 */

#ifndef BOINK_INTERFACE_HH
#define BOINK_INTERFACE_HH

// make sure that <condition_variable> does
// not included from the standard lib because
// it breaks conda builds with cppyy=1.5.5
#define _GLIBCXX_CONDITION_VARIABLE 1

#include "boink/boink.hh"

#include "boink/storage/storage.hh"
#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/partitioned_storage.hh"
#include "boink/storage/sparseppstorage.hh"
#include "boink/storage/sparsepp/spp.h"

#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/hashshifter.hh"
#include "boink/hashing/hashextender.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/ukhshashshifter.hh"

#include "boink/kmers/kmerclient.hh"
#include "boink/sequences/alphabets.hh"
#include "boink/sequences/exceptions.hh"

#include "boink/parsing/parsing.hh"
#include "boink/parsing/readers.hh"

//#include "boink/events.hh"
//#include "boink/event_types.hh"

#include "boink/metrics.hh"


#include "boink/dbg.hh"
#include "boink/pdbg.hh"
#include "boink/traversal.hh"

#include "boink/reporting/reporters.hh"

#include "boink/processors.hh"

#include "boink/cdbg/cdbg_types.hh"
#include "boink/cdbg/compactor.hh"
#include "boink/cdbg/cdbg.hh"
#include "boink/cdbg/metrics.hh"
#include "boink/cdbg/udbg.hh"
#include "boink/cdbg/utagger.hh"
#include "boink/cdbg/saturating_compactor.hh"
//#include "boink/cdbg/ucompactor.hh"

#include "boink/minimizers.hh"
#include "boink/signatures/ukhs_signature.hh"
#include "boink/signatures/sourmash_signature.hh"

#include "boink/benchmarks/bench_storage.hh"

#include <set>

#define STLTYPES_EXPLICIT_INSTANTIATION_DECL(STLTYPE, TTYPE)                      \
template class std::STLTYPE< TTYPE >;                                             \
template class __gnu_cxx::__normal_iterator<TTYPE*, std::STLTYPE< TTYPE > >;      \
template class __gnu_cxx::__normal_iterator<const TTYPE*, std::STLTYPE< TTYPE > >;\
namespace __gnu_cxx {                                                             \
template bool operator==(const std::STLTYPE< TTYPE >::iterator&,                  \
                         const std::STLTYPE< TTYPE >::iterator&);                 \
template bool operator!=(const std::STLTYPE< TTYPE >::iterator&,                  \
                         const std::STLTYPE< TTYPE >::iterator&);                 \
}

namespace boink {

//
// Readers
//

extern template class parsing::FastxParser<DNA_SIMPLE>;
extern template class parsing::SequenceReader<parsing::FastxParser<DNA_SIMPLE>>;
extern template class parsing::SplitPairedReader<parsing::FastxParser<DNA_SIMPLE>>;

namespace hashing {

    //
    // Kmer Models
    //

    extern template class HashModel<uint64_t>;
    extern template class CanonicalModel<HashModel<uint64_t>>;

    extern template class Partitioned<HashModel<uint64_t>>;
    extern template class Partitioned<CanonicalModel<uint64_t>>;

    extern template class KmerModel<HashModel<uint64_t>>;
    extern template class KmerModel<CanonicalModel<HashModel<uint64_t>>>;
    extern template class KmerModel<UnikmerWmer>;
    extern template class KmerModel<CanonicalUnikmerWmer>;

    extern template class ShiftModel<Hash, DIR_LEFT>;
    extern template class ShiftModel<Hash, DIR_RIGHT>;
    extern template class ShiftModel<Canonical, DIR_LEFT>;
    extern template class ShiftModel<Canonical, DIR_RIGHT>;

    //
    // Shifters
    //

    extern template class RollingHashShifter<HashModel<uint64_t>>;
    extern template class RollingHashShifter<CanonicalModel<uint64_t>>;

    extern template class UnikmerShifter<FwdRollingShifter>;
    extern template class UnikmerShifter<CanRollingShifter>;

    //
    // HashExtenders
    //

    extern template class HashExtender<FwdRollingShifter>;
    extern template class HashExtender<CanRollingShifter>;

    extern template class HashExtender<FwdUnikmerShifter>;
    extern template class HashExtender<CanUnikmerShifter>;

    //
    // KmerIterator
    //

    extern template class KmerIterator<FwdRollingShifter>;
    extern template class KmerIterator<CanRollingShifter>;

    extern template class KmerIterator<FwdUnikmerShifter>;
    extern template class KmerIterator<CanUnikmerShifter>;

    extern template class KmerIterator<HashExtender<FwdRollingShifter>>;
    extern template class KmerIterator<HashExtender<CanRollingShifter>>;

    extern template class KmerIterator<HashExtender<FwdUnikmerShifter>>;
    extern template class KmerIterator<HashExtender<CanUnikmerShifter>>;

    extern template class KmerSpanMixinImpl<false>;
    extern template class KmerSpanMixinImpl<true>;

    //
    // UKHS
    //

    extern template class UKHS<FwdRollingShifter>;
    extern template class UKHS<CanRollingShifter>;


}

//
// dBG
//

extern template class dBG<storage::BitStorage, hashing::FwdRollingShifter>;
extern template class dBG<storage::BitStorage, hashing::CanRollingShifter>;
extern template class dBG<storage::BitStorage, hashing::FwdUnikmerShifter>;
extern template class dBG<storage::BitStorage, hashing::CanUnikmerShifter>;

extern template class dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>;
extern template class dBG<storage::SparseppSetStorage, hashing::CanRollingShifter>;
extern template class dBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>;
extern template class dBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>;

extern template class dBG<storage::ByteStorage, hashing::FwdRollingShifter>;
extern template class dBG<storage::ByteStorage, hashing::CanRollingShifter>;
extern template class dBG<storage::ByteStorage, hashing::FwdUnikmerShifter>;
extern template class dBG<storage::ByteStorage, hashing::CanUnikmerShifter>;

extern template class dBG<storage::NibbleStorage, hashing::FwdRollingShifter>;
extern template class dBG<storage::NibbleStorage, hashing::CanRollingShifter>;
extern template class dBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>;
extern template class dBG<storage::NibbleStorage, hashing::CanUnikmerShifter>;

extern template class dBG<storage::QFStorage, hashing::FwdRollingShifter>;
extern template class dBG<storage::QFStorage, hashing::CanRollingShifter>;
extern template class dBG<storage::QFStorage, hashing::FwdUnikmerShifter>;
extern template class dBG<storage::QFStorage, hashing::CanUnikmerShifter>;

//
//  PdBG
//

extern template class PdBG<storage::BitStorage, hashing::FwdRollingShifter>;
extern template class PdBG<storage::BitStorage, hashing::CanRollingShifter>;

extern template class PdBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>;
extern template class PdBG<storage::SparseppSetStorage, hashing::CanRollingShifter>;

extern template class PdBG<storage::ByteStorage, hashing::FwdRollingShifter>;
extern template class PdBG<storage::ByteStorage, hashing::CanRollingShifter>;

extern template class PdBG<storage::NibbleStorage, hashing::FwdRollingShifter>;
extern template class PdBG<storage::NibbleStorage, hashing::CanRollingShifter>;

extern template class PdBG<storage::QFStorage, hashing::FwdRollingShifter>;
extern template class PdBG<storage::QFStorage, hashing::CanRollingShifter>;

//
// Walkers
//

extern template class dBGWalker<dBG<storage::BitStorage, hashing::FwdRollingShifter>>;
extern template class dBGWalker<dBG<storage::BitStorage, hashing::CanRollingShifter>>;
extern template class dBGWalker<dBG<storage::BitStorage, hashing::FwdUnikmerShifter>>;
extern template class dBGWalker<dBG<storage::BitStorage, hashing::CanUnikmerShifter>>;

extern template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>>;
extern template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::CanRollingShifter>>;
extern template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>>;
extern template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>>;

extern template class dBGWalker<dBG<storage::ByteStorage, hashing::FwdRollingShifter>>;
extern template class dBGWalker<dBG<storage::ByteStorage, hashing::CanRollingShifter>>;
extern template class dBGWalker<dBG<storage::ByteStorage, hashing::FwdUnikmerShifter>>;
extern template class dBGWalker<dBG<storage::ByteStorage, hashing::CanUnikmerShifter>>;

extern template class dBGWalker<dBG<storage::NibbleStorage, hashing::FwdRollingShifter>>;
extern template class dBGWalker<dBG<storage::NibbleStorage, hashing::CanRollingShifter>>;
extern template class dBGWalker<dBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>>;
extern template class dBGWalker<dBG<storage::NibbleStorage, hashing::CanUnikmerShifter>>;

extern template class dBGWalker<dBG<storage::QFStorage, hashing::FwdRollingShifter>>;
extern template class dBGWalker<dBG<storage::QFStorage, hashing::CanRollingShifter>>;
extern template class dBGWalker<dBG<storage::QFStorage, hashing::FwdUnikmerShifter>>;
extern template class dBGWalker<dBG<storage::QFStorage, hashing::CanUnikmerShifter>>;

extern template class hashing::KmerIterator<dBGWalker<dBG<storage::BitStorage, hashing::FwdRollingShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::BitStorage, hashing::CanRollingShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::BitStorage, hashing::FwdUnikmerShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::BitStorage, hashing::CanUnikmerShifter>>>;

extern template class hashing::KmerIterator<dBGWalker<dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::SparseppSetStorage, hashing::CanRollingShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>>>;

extern template class hashing::KmerIterator<dBGWalker<dBG<storage::ByteStorage, hashing::FwdRollingShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::ByteStorage, hashing::CanRollingShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::ByteStorage, hashing::FwdUnikmerShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::ByteStorage, hashing::CanUnikmerShifter>>>;

extern template class hashing::KmerIterator<dBGWalker<dBG<storage::NibbleStorage, hashing::FwdRollingShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::NibbleStorage, hashing::CanRollingShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::NibbleStorage, hashing::CanUnikmerShifter>>>;

extern template class hashing::KmerIterator<dBGWalker<dBG<storage::QFStorage, hashing::FwdRollingShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::QFStorage, hashing::CanRollingShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::QFStorage, hashing::FwdUnikmerShifter>>>;
extern template class hashing::KmerIterator<dBGWalker<dBG<storage::QFStorage, hashing::CanUnikmerShifter>>>;


//
// cDBG
//

extern template class cdbg::cDBG<dBG<storage::BitStorage, hashing::FwdRollingShifter>>;
extern template class cdbg::cDBG<dBG<storage::BitStorage, hashing::CanRollingShifter>>;

extern template class cdbg::cDBG<dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>>;
extern template class cdbg::cDBG<dBG<storage::SparseppSetStorage, hashing::CanRollingShifter>>;

extern template class cdbg::cDBG<dBG<storage::ByteStorage, hashing::FwdRollingShifter>>;
extern template class cdbg::cDBG<dBG<storage::ByteStorage, hashing::CanRollingShifter>>;

extern template class cdbg::cDBG<dBG<storage::NibbleStorage, hashing::FwdRollingShifter>>;
extern template class cdbg::cDBG<dBG<storage::NibbleStorage, hashing::CanRollingShifter>>;

extern template class cdbg::cDBG<dBG<storage::QFStorage, hashing::FwdRollingShifter>>;
extern template class cdbg::cDBG<dBG<storage::QFStorage, hashing::CanRollingShifter>>;

//
// uDBG
//

extern template class cdbg::uDBG<storage::BitStorage>;
extern template class cdbg::uDBG<storage::ByteStorage>;
extern template class cdbg::uDBG<storage::NibbleStorage>;
extern template class cdbg::uDBG<storage::QFStorage>;
extern template class cdbg::uDBG<storage::SparseppSetStorage>;

//
// UTagger
//

extern template class cdbg::UTagger<storage::BitStorage>;
extern template class cdbg::UTagger<storage::ByteStorage>;
extern template class cdbg::UTagger<storage::NibbleStorage>;
extern template class cdbg::UTagger<storage::QFStorage>;
extern template class cdbg::UTagger<storage::SparseppSetStorage>;

//
// StreamingCompactor
//

extern template class cdbg::StreamingCompactor<dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>>;
extern template class cdbg::StreamingCompactor<dBG<storage::BitStorage, hashing::FwdRollingShifter>>;
extern template class cdbg::StreamingCompactor<dBG<storage::ByteStorage, hashing::FwdRollingShifter>>;
extern template class cdbg::StreamingCompactor<dBG<storage::NibbleStorage, hashing::FwdRollingShifter>>;
extern template class cdbg::StreamingCompactor<dBG<storage::QFStorage, hashing::FwdRollingShifter>>;

//
// UnikmerSignature
//

extern template class signatures::UnikmerSignature<storage::BitStorage, hashing::FwdRollingShifter>;
extern template class signatures::UnikmerSignature<storage::BitStorage, hashing::CanRollingShifter>;

extern template class signatures::UnikmerSignature<storage::SparseppSetStorage, hashing::FwdRollingShifter>;
extern template class signatures::UnikmerSignature<storage::SparseppSetStorage, hashing::CanRollingShifter>;

extern template class signatures::UnikmerSignature<storage::ByteStorage, hashing::FwdRollingShifter>;
extern template class signatures::UnikmerSignature<storage::ByteStorage, hashing::CanRollingShifter>;

extern template class signatures::UnikmerSignature<storage::NibbleStorage, hashing::FwdRollingShifter>;
extern template class signatures::UnikmerSignature<storage::NibbleStorage, hashing::CanRollingShifter>;

extern template class signatures::UnikmerSignature<storage::QFStorage, hashing::FwdRollingShifter>;
extern template class signatures::UnikmerSignature<storage::QFStorage, hashing::CanRollingShifter>;

//
// Minimizers
//

extern template class InteriorMinimizer<uint64_t>;
extern template class WKMinimizer<hashing::FwdRollingShifter>;
extern template class WKMinimizer<hashing::CanRollingShifter>;


}

#endif
