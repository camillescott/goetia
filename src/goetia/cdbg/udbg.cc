/* goetia/cdbg/udbg.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "goetia/cdbg/udbg.hh"
#include "goetia/storage/storage_types.hh"

namespace goetia {

    template class uDBG<BitStorage>;
    template class uDBG<ByteStorage>;
    template class uDBG<NibbleStorage>;
    template class uDBG<QFStorage>;
    template class uDBG<SparseppSetStorage>;

}
