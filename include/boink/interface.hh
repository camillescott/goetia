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
#include "boink/hashing/unikmershifter.hh"
#include "boink/hashing/canonical.hh"

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
#include "boink/cdbg/usparsetagger.hh"
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

#endif
