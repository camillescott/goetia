/* boink.hh
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "boink/processors.hh"

#include "boink/dbg.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/sparseppstorage.hh"

#include "sourmash/kmer_min_hash.hh"

namespace boink {

SourmashSignatureProcessor
::SourmashSignatureProcessor(KmerMinHash * signature,
                             uint64_t fine_interval,
                             uint64_t medium_interval,
                             uint64_t coarse_interval)
    : Base(fine_interval, medium_interval, coarse_interval),
      signature(signature)
{
}

void SourmashSignatureProcessor::process_sequence(const parsing::Read& read) {
    signature->add_sequence(read.cleaned_seq.c_str(), false);
}

void SourmashSignatureProcessor::report() {

}

}

template class boink::FileConsumer<boink::dBG<boink::storage::BitStorage,
                                              boink::hashing::RollingHashShifter>>;
template class boink::FileConsumer<boink::dBG<boink::storage::ByteStorage,
                                              boink::hashing::RollingHashShifter>>;
template class boink::FileConsumer<boink::dBG<boink::storage::NibbleStorage,
                                              boink::hashing::RollingHashShifter>>;
template class boink::FileConsumer<boink::dBG<boink::storage::QFStorage,
                                              boink::hashing::RollingHashShifter>>;
template class boink::FileConsumer<boink::dBG<boink::storage::SparseppSetStorage,
                                              boink::hashing::RollingHashShifter>>;


template class boink::DecisionNodeProcessor<boink::dBG<boink::storage::BitStorage,
                                                       boink::hashing::RollingHashShifter>>;
template class boink::DecisionNodeProcessor<boink::dBG<boink::storage::ByteStorage,
                                                       boink::hashing::RollingHashShifter>>;
template class boink::DecisionNodeProcessor<boink::dBG<boink::storage::NibbleStorage,
                                                       boink::hashing::RollingHashShifter>>;
template class boink::DecisionNodeProcessor<boink::dBG<boink::storage::QFStorage,
                                                       boink::hashing::RollingHashShifter>>;
template class boink::DecisionNodeProcessor<boink::dBG<boink::storage::SparseppSetStorage,
                                                       boink::hashing::RollingHashShifter>>;


template class boink::StreamingCompactorProcessor<boink::dBG<boink::storage::BitStorage,
                                                             boink::hashing::RollingHashShifter>>;
template class boink::StreamingCompactorProcessor<boink::dBG<boink::storage::ByteStorage,
                                                             boink::hashing::RollingHashShifter>>;
template class boink::StreamingCompactorProcessor<boink::dBG<boink::storage::NibbleStorage,
                                                             boink::hashing::RollingHashShifter>>;
template class boink::StreamingCompactorProcessor<boink::dBG<boink::storage::QFStorage,
                                                             boink::hashing::RollingHashShifter>>;
template class boink::StreamingCompactorProcessor<boink::dBG<boink::storage::SparseppSetStorage,
                                                             boink::hashing::RollingHashShifter>>;

