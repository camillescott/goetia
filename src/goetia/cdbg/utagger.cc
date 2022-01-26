/**
 * (c) Camille Scott, 2019
 * File   : utagger.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 01.11.2019
 */

#include "goetia/cdbg/utagger.hh"
#include "goetia/storage/storage_types.hh"

namespace goetia {

    template class UTagger<BitStorage>;
    template class UTagger<ByteStorage>;
    template class UTagger<NibbleStorage>;
    template class UTagger<QFStorage>;
    template class UTagger<SparseppSetStorage>;

}
