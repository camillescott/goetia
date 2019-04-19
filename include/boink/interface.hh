#include "boink/boink.hh"

#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/hashshifter.hh"
#include "boink/hashing/exceptions.hh"
#include "boink/hashing/alphabets.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/hashing/rollinghashshifter.hh"

#include "boink/kmers/kmerclient.hh"

#include "boink/parsing/parsing.hh"
#include "boink/parsing/readers.hh"

#include "boink/events.hh"
#include "boink/event_types.hh"

#include "boink/metrics.hh"

#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/partitioned_storage.hh"
#include "boink/storage/sparseppstorage.hh"

#include "boink/dbg.hh"
#include "boink/pdbg.hh"
#include "boink/assembly.hh"

#include "boink/reporting/cdbg_writer_reporter.hh"
#include "boink/reporting/report_types.hh"
#include "boink/reporting/reporters.hh"
#include "boink/reporting/cdbg_component_reporter.hh"
#include "boink/reporting/cdbg_history_reporter.hh"
#include "boink/reporting/streaming_compactor_reporter.hh"
#include "boink/reporting/cdbg_unitig_reporter.hh"

#include "boink/processors.hh"

#include "boink/normalization/diginorm.hh"

#include "boink/cdbg/cdbg_types.hh"
#include "boink/cdbg/compactor.hh"
#include "boink/cdbg/cdbg.hh"
#include "boink/cdbg/solid_compactor.hh"
#include "boink/cdbg/metrics.hh"

#include "boink/minimizers.hh"
#include "boink/signatures/ukhs_signature.hh"

#include "prometheus/exposer.h"

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

template <typename ReaderType=boink::parsing::FastxReader>
struct ParserFactory {
	static std::shared_ptr<boink::parsing::ReadParser<boink::parsing::FastxReader>> build(const std::string& filename) {
		return boink::parsing::get_parser<boink::parsing::FastxReader>(filename);
	}
};

template struct ParserFactory<boink::parsing::FastxReader>;

/*
 * BitStorage  declarations
 */

// dBG
template class boink::dBG<boink::storage::BitStorage,
                          boink::hashing::RollingHashShifter>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::BitStorage,
                                            boink::hashing::RollingHashShifter>>;

// PdBG
template class boink::PdBG<boink::storage::BitStorage>;
template class std::enable_shared_from_this<boink::PdBG<boink::storage::BitStorage>>;

// Assemblers, Compactors
template class boink::AssemblerMixin<boink::dBG<boink::storage::BitStorage,
                                                 boink::hashing::RollingHashShifter>>;
template class boink::CompactorMixin<boink::dBG<boink::storage::BitStorage,
                                                 boink::hashing::RollingHashShifter>>;
template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::BitStorage,
                                               boink::hashing::RollingHashShifter>>;

//Processors
template class boink::FileConsumer<boink::dBG<boink::storage::BitStorage,
                                              boink::hashing::RollingHashShifter>>;
template class boink::DecisionNodeProcessor<boink::dBG<boink::storage::BitStorage,
                                                       boink::hashing::RollingHashShifter>>;
template class boink::StreamingCompactorProcessor<boink::dBG<boink::storage::BitStorage,
                                                             boink::hashing::RollingHashShifter>>;

//Signatures
extern template class boink::signatures::UnikmerSignature<boink::storage::BitStorage>;

/*
 * ByteStorage  declarations
 */

// dBG
template class boink::dBG<boink::storage::ByteStorage,
                          boink::hashing::RollingHashShifter>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::ByteStorage,
                                            boink::hashing::RollingHashShifter>>;

// PdBG
template class boink::PdBG<boink::storage::ByteStorage>;
template class std::enable_shared_from_this<boink::PdBG<boink::storage::ByteStorage>>;

// Assemblers, Compactors
template class boink::AssemblerMixin<boink::dBG<boink::storage::ByteStorage,
                                                 boink::hashing::RollingHashShifter>>;
template class boink::CompactorMixin<boink::dBG<boink::storage::ByteStorage,
                                                 boink::hashing::RollingHashShifter>>;
template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::ByteStorage,
                                               boink::hashing::RollingHashShifter>>;

//Processors
template class boink::FileConsumer<boink::dBG<boink::storage::ByteStorage,
                                              boink::hashing::RollingHashShifter>>;
template class boink::DecisionNodeProcessor<boink::dBG<boink::storage::ByteStorage,
                                                       boink::hashing::RollingHashShifter>>;
template class boink::StreamingCompactorProcessor<boink::dBG<boink::storage::ByteStorage,
                                                             boink::hashing::RollingHashShifter>>;

//Signatures
extern template class boink::signatures::UnikmerSignature<boink::storage::ByteStorage>;

/*
 * NibbleStorage  declarations
 */

// dBG
template class boink::dBG<boink::storage::NibbleStorage,
                          boink::hashing::RollingHashShifter>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::NibbleStorage,
                                            boink::hashing::RollingHashShifter>>;

// PdBG
template class boink::PdBG<boink::storage::NibbleStorage>;
template class std::enable_shared_from_this<boink::PdBG<boink::storage::NibbleStorage>>;

// Assemblers, Compactors
template class boink::AssemblerMixin<boink::dBG<boink::storage::NibbleStorage,
                                                 boink::hashing::RollingHashShifter>>;
template class boink::CompactorMixin<boink::dBG<boink::storage::NibbleStorage,
                                                 boink::hashing::RollingHashShifter>>;
template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::NibbleStorage,
                                               boink::hashing::RollingHashShifter>>;

//Processors
template class boink::FileConsumer<boink::dBG<boink::storage::NibbleStorage,
                                              boink::hashing::RollingHashShifter>>;
template class boink::DecisionNodeProcessor<boink::dBG<boink::storage::NibbleStorage,
                                                       boink::hashing::RollingHashShifter>>;
template class boink::StreamingCompactorProcessor<boink::dBG<boink::storage::NibbleStorage,
                                                             boink::hashing::RollingHashShifter>>;

//Signatures
extern template class boink::signatures::UnikmerSignature<boink::storage::NibbleStorage>;

/*
 * QFStorage  declarations
 */

// dBG
template class boink::dBG<boink::storage::QFStorage,
                          boink::hashing::RollingHashShifter>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::QFStorage,
                                            boink::hashing::RollingHashShifter>>;

// PdBG
template class boink::PdBG<boink::storage::QFStorage>;
template class std::enable_shared_from_this<boink::PdBG<boink::storage::QFStorage>>;

// Assemblers, Compactors
template class boink::AssemblerMixin<boink::dBG<boink::storage::QFStorage,
                                                 boink::hashing::RollingHashShifter>>;
template class boink::CompactorMixin<boink::dBG<boink::storage::QFStorage,
                                                 boink::hashing::RollingHashShifter>>;
template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::QFStorage,
                                               boink::hashing::RollingHashShifter>>;

//Processors
template class boink::FileConsumer<boink::dBG<boink::storage::QFStorage,
                                              boink::hashing::RollingHashShifter>>;
template class boink::DecisionNodeProcessor<boink::dBG<boink::storage::QFStorage,
                                                       boink::hashing::RollingHashShifter>>;
template class boink::StreamingCompactorProcessor<boink::dBG<boink::storage::QFStorage,
                                                             boink::hashing::RollingHashShifter>>;

//Signatures
extern template class boink::signatures::UnikmerSignature<boink::storage::QFStorage>;

/*
 * SparseppSetStorage  declarations
 */

// dBG
template class boink::dBG<boink::storage::SparseppSetStorage,
                          boink::hashing::RollingHashShifter>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::SparseppSetStorage,
                                            boink::hashing::RollingHashShifter>>;

// PdBG
template class boink::PdBG<boink::storage::SparseppSetStorage>;
template class std::enable_shared_from_this<boink::PdBG<boink::storage::SparseppSetStorage>>;

// Assemblers, Compactors
template class boink::AssemblerMixin<boink::dBG<boink::storage::SparseppSetStorage,
                                                 boink::hashing::RollingHashShifter>>;
template class boink::CompactorMixin<boink::dBG<boink::storage::SparseppSetStorage,
                                                 boink::hashing::RollingHashShifter>>;
template class boink::cdbg::StreamingCompactor<boink::dBG<boink::storage::SparseppSetStorage,
                                               boink::hashing::RollingHashShifter>>;

//Processors
template class boink::FileConsumer<boink::dBG<boink::storage::SparseppSetStorage,
                                              boink::hashing::RollingHashShifter>>;
template class boink::DecisionNodeProcessor<boink::dBG<boink::storage::SparseppSetStorage,
                                                       boink::hashing::RollingHashShifter>>;
template class boink::StreamingCompactorProcessor<boink::dBG<boink::storage::SparseppSetStorage,
                                                             boink::hashing::RollingHashShifter>>;

//Signatures
extern template class boink::signatures::UnikmerSignature<boink::storage::SparseppSetStorage>;
