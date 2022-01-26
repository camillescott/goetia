/* metrics.cc -- utilities for metrics gathering
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "goetia/metrics.hh"

namespace goetia {

    template class ReservoirSample<uint64_t>;
    template class ReservoirSample<double>;
    template class ReservoirSample<float>;
    template class ReservoirSample<uint32_t>;

}
