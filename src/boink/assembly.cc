#include "boink/assembly.hh"

#include "boink/dbg.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/sparseppstorage.hh"

using namespace boink;


template class boink::AssemblerMixin<boink::dBG<boink::storage::BitStorage,
                                                boink::hashing::RollingHashShifter>>;
template class boink::AssemblerMixin<boink::dBG<boink::storage::ByteStorage,
                                                boink::hashing::RollingHashShifter>>;
template class boink::AssemblerMixin<boink::dBG<boink::storage::NibbleStorage,
                                                boink::hashing::RollingHashShifter>>;
template class boink::AssemblerMixin<boink::dBG<boink::storage::QFStorage,
                                                boink::hashing::RollingHashShifter>>;
template class boink::AssemblerMixin<boink::dBG<boink::storage::SparseppSetStorage,
                                                boink::hashing::RollingHashShifter>>;


template class boink::CompactorMixin<boink::dBG<boink::storage::BitStorage,
                                                boink::hashing::RollingHashShifter>>;
template class boink::CompactorMixin<boink::dBG<boink::storage::ByteStorage,
                                                boink::hashing::RollingHashShifter>>;
template class boink::CompactorMixin<boink::dBG<boink::storage::NibbleStorage,
                                                boink::hashing::RollingHashShifter>>;
template class boink::CompactorMixin<boink::dBG<boink::storage::QFStorage,
                                                boink::hashing::RollingHashShifter>>;
template class boink::CompactorMixin<boink::dBG<boink::storage::SparseppSetStorage,
                                                boink::hashing::RollingHashShifter>>;
