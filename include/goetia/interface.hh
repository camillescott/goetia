/**
 * (c) Camille Scott, 2019
 * File   : interface.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 25.09.2019
 */

#ifndef GOETIA_INTERFACE_HH
#define GOETIA_INTERFACE_HH

// make sure that <condition_variable> does
// not included from the standard lib because
// it breaks conda builds with cppyy=1.5.5
#define _GLIBCXX_CONDITION_VARIABLE 1

#include "goetia/goetia.hh"

#include "goetia/storage/storage.hh"
#include "goetia/storage/nibblestorage.hh"
#include "goetia/storage/bitstorage.hh"
#include "goetia/storage/qfstorage.hh"
#include "goetia/storage/bytestorage.hh"
#include "goetia/storage/partitioned_storage.hh"
#include "goetia/storage/sparseppstorage.hh"
#include "goetia/storage/phmapstorage.hh"
#include "goetia/storage/sparsepp/spp.h"

#include "goetia/hashing/kmeriterator.hh"
#include "goetia/hashing/hashshifter.hh"
#include "goetia/hashing/hashextender.hh"
#include "goetia/hashing/ukhs.hh"
#include "goetia/hashing/rollinghashshifter.hh"
#include "goetia/hashing/unikmershifter.hh"
#include "goetia/hashing/canonical.hh"

#include "goetia/sequences/alphabets.hh"
#include "goetia/sequences/exceptions.hh"

#include "goetia/parsing/parsing.hh"
#include "goetia/parsing/readers.hh"

//#include "goetia/events.hh"
//#include "goetia/event_types.hh"

#include "goetia/metrics.hh"


#include "goetia/dbg.hh"
#include "goetia/pdbg.hh"
#include "goetia/traversal.hh"
#include "goetia/solidifier.hh"

#include "goetia/processors.hh"

#include "goetia/cdbg/cdbg_types.hh"
#include "goetia/cdbg/compactor.hh"
#include "goetia/cdbg/cdbg.hh"
#include "goetia/cdbg/metrics.hh"
#include "goetia/cdbg/udbg.hh"
#include "goetia/cdbg/utagger.hh"
#include "goetia/cdbg/saturating_compactor.hh"
//#include "goetia/cdbg/ucompactor.hh"

#include "goetia/minimizers.hh"
#include "goetia/signatures/ukhs_signature.hh"
#include "goetia/signatures/sourmash_signature.hh"
#include "goetia/signatures/sourmash/sourmash.hpp"

#include "goetia/benchmarks/bench_storage.hh"

//#include <set>

//#define STLTYPES_EXPLICIT_INSTANTIATION_DECL(STLTYPE, TTYPE)                      \
//template class std::STLTYPE< TTYPE >;                                             \
//template class __gnu_cxx::__normal_iterator<TTYPE*, std::STLTYPE< TTYPE > >;      \
//template class __gnu_cxx::__normal_iterator<const TTYPE*, std::STLTYPE< TTYPE > >;\
//namespace __gnu_cxx {                                                             \
//template bool operator==(const std::STLTYPE< TTYPE >::iterator&,                  \
//                         const std::STLTYPE< TTYPE >::iterator&);                 \
//template bool operator!=(const std::STLTYPE< TTYPE >::iterator&,                  \
//                         const std::STLTYPE< TTYPE >::iterator&);                 \
//}

#endif
