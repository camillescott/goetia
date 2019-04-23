#include "boink/pdbg.hh"

#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/sparseppstorage.hh"

#include <memory>

using namespace boink;


template class boink::PdBG<boink::storage::BitStorage>;
template class boink::PdBG<boink::storage::ByteStorage>;
template class boink::PdBG<boink::storage::NibbleStorage>;
template class boink::PdBG<boink::storage::QFStorage>;
template class boink::PdBG<boink::storage::SparseppSetStorage>;

