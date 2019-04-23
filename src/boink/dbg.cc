#include "boink/dbg.hh"

#include "boink/hashing/rollinghashshifter.hh"
#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/sparseppstorage.hh"

#include <memory>

using namespace boink;

template class boink::dBG<boink::storage::BitStorage,
                          boink::hashing::RollingHashShifter>;
template class boink::dBG<boink::storage::ByteStorage,
                          boink::hashing::RollingHashShifter>;
template class boink::dBG<boink::storage::NibbleStorage,
                          boink::hashing::RollingHashShifter>;
template class boink::dBG<boink::storage::QFStorage,
                          boink::hashing::RollingHashShifter>;
template class boink::dBG<boink::storage::SparseppSetStorage,
                          boink::hashing::RollingHashShifter>;

template class std::enable_shared_from_this<boink::dBG<boink::storage::BitStorage,
                                                       boink::hashing::RollingHashShifter>>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::ByteStorage,
                                                       boink::hashing::RollingHashShifter>>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::NibbleStorage,
                                                       boink::hashing::RollingHashShifter>>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::QFStorage,
                                                       boink::hashing::RollingHashShifter>>;
template class std::enable_shared_from_this<boink::dBG<boink::storage::SparseppSetStorage,
                                                       boink::hashing::RollingHashShifter>>;


