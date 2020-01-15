/* boink.hh
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "boink/minimizers.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/ukhshashshifter.hh"

namespace boink {

template class WKMinimizer<hashing::FwdRollingShifter>;
template class WKMinimizer<hashing::CanRollingShifter>;

}
