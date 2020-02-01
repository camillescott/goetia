/* rollinghashshifter.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "boink/hashing/rollinghashshifter.hh"

namespace boink {

namespace hashing {

template class RollingHashShifter<Hash<uint64_t>>;
template class RollingHashShifter<Canonical<uint64_t>>;

}
}

