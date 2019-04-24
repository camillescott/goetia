/* kmerclient.hh: base class for things parametrized by K
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_KMERCLIENT_HH
#define BOINK_KMERCLIENT_HH

#include <cstdint>

namespace boink {
namespace kmers {

class KmerClient {
protected:
    const uint16_t _K;
    explicit KmerClient(uint16_t K) : _K(K) {}

public:
    const uint16_t K() const { return _K; }
};

}
}

#endif
