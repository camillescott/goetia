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

    template class cdbg::UTagger<storage::BitStorage>;
    template class cdbg::UTagger<storage::ByteStorage>;
    template class cdbg::UTagger<storage::NibbleStorage>;
    template class cdbg::UTagger<storage::QFStorage>;
    template class cdbg::UTagger<storage::SparseppSetStorage>;

}
