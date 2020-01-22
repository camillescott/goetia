/* boink/cdbg/udbg.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "boink/cdbg/udbg.hh"
#include "boink/storage/storage_types.hh"

namespace boink {

    template class cdbg::uDBG<storage::BitStorage>;
    template class cdbg::uDBG<storage::ByteStorage>;
    template class cdbg::uDBG<storage::NibbleStorage>;
    template class cdbg::uDBG<storage::QFStorage>;
    template class cdbg::uDBG<storage::SparseppSetStorage>;

}
