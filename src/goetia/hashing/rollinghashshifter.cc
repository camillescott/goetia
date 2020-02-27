/* rollinghashshifter.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "goetia/hashing/rollinghashshifter.hh"

namespace goetia {

namespace hashing {

template class LemireShifterPolicy<Hash<uint64_t>>;
template class LemireShifterPolicy<Canonical<uint64_t>>;

}
}

