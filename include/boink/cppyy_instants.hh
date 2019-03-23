#include "boink/hashing/alphabets.hh"
#include "boink/hashing/exceptions.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/hashshifter.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/rollinghashshifter.hh"

#include "boink/assembly.hh"

#include "boink/dbg.hh"

#include "boink/storage/storage.hh"
#include "boink/storage/sparseppstorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/nibblestorage.hh"
#include "boink/storage/qfstorage.hh"

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

//#define STLTYPES_EXPLICIT_INSTANTIATION_DECL(STLTYPE, TTYPE)                        \
template class std::STLTYPE< TTYPE >;                                               \
template class __gnu_cxx::__normal_iterator<TTYPE*, std::STLTYPE< TTYPE > >;\
template class __gnu_cxx::__normal_iterator<const TTYPE*, std::STLTYPE< TTYPE > >; \
namespace __gnu_cxx {                                                           \
template bool operator==(const __normal_iterator<TTYPE*, std::STLTYPE< TTYPE > >&,  \
                         const __normal_iterator<TTYPE*, std::STLTYPE< TTYPE > >&); \
template bool operator!=(const __normal_iterator<TTYPE*, std::STLTYPE< TTYPE > >&,  \
                         const __normal_iterator<TTYPE*, std::STLTYPE< TTYPE > >&); \
}
//template bool operator!=(const std::STLTYPE< TTYPE >::iterator&,                  \
                         const std::STLTYPE< TTYPE >::iterator&);                 \




STLTYPES_EXPLICIT_INSTANTIATION_DECL(vector, boink::hashing::hash_t)
STLTYPES_EXPLICIT_INSTANTIATION_DECL(vector, boink::hashing::shift_t)
STLTYPES_EXPLICIT_INSTANTIATION_DECL(vector, boink::hashing::kmer_t)


template class std::pair<boink::hashing::shift_t, boink::hashing::shift_t>;

template class boink::hashing::KmerIterator<boink::hashing::RollingHashShifter>;

//
// dBG storage specializations
//


// SparseppSetStorage
template class boink::dBG<boink::storage::SparseppSetStorage, 
                          boink::hashing::RollingHashShifter>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::SparseppSetStorage, 
                                                       boink::hashing::RollingHashShifter>>;
template class boink::AssemblerMixin<boink::dBG<boink::storage::SparseppSetStorage,
                                                boink::hashing::RollingHashShifter>>;

// BitStorage
template class boink::dBG<boink::storage::BitStorage,
                          boink::hashing::RollingHashShifter>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::BitStorage,
                                            boink::hashing::RollingHashShifter>>;
template class boink::AssemblerMixin<boink::dBG<boink::storage::BitStorage,
                                                 boink::hashing::RollingHashShifter>>;

// ByteStorage
template class boink::dBG<boink::storage::ByteStorage,
                          boink::hashing::RollingHashShifter>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::ByteStorage,
                                            boink::hashing::RollingHashShifter>>;
template class boink::AssemblerMixin<boink::dBG<boink::storage::ByteStorage,
                                                 boink::hashing::RollingHashShifter>>;

// NibbleStorage
template class boink::dBG<boink::storage::NibbleStorage,
                          boink::hashing::RollingHashShifter>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::NibbleStorage,
                                            boink::hashing::RollingHashShifter>>;
template class boink::AssemblerMixin<boink::dBG<boink::storage::NibbleStorage,
                                                 boink::hashing::RollingHashShifter>>;

// QFStorage
template class boink::dBG<boink::storage::QFStorage,
                          boink::hashing::RollingHashShifter>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::QFStorage,
                                            boink::hashing::RollingHashShifter>>;
template class boink::AssemblerMixin<boink::dBG<boink::storage::QFStorage,
                                                 boink::hashing::RollingHashShifter>>;
