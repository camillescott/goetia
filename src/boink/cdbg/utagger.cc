/**
 * (c) Camille Scott, 2019
 * File   : utagger.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 01.11.2019
 */

#include "boink/cdbg/utagger.hh"
#include "boink/storage/storage_types.hh"

namespace boink {
namespace cdbg {

}
}

template class boink::cdbg::UTagger<boink::storage::BitStorage>;
//template class boink::cdbg::UTagger<boink::storage::ByteStorage>;
//template class boink::cdbg::UTagger<boink::storage::NibbleStorage>;
//template class boink::cdbg::UTagger<boink::storage::QFStorage>;
//template class boink::cdbg::UTagger<boink::storage::SparseppSetStorage>;
