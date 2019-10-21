/* metrics.hh -- utilities for metrics gathering
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_UTILITY_METRICS_HH
#define BOINK_UTILITY_METRICS_HH

#include <algorithm>
#include <atomic>
#include <iterator>
#include <random>

namespace boink {
namespace metrics {

/**
 * @Synopsis  Streaming reservoir sampling: continuously calling
 *            sample(item) on a stream of N items will eventually converge
 *            to sampling each item with probability sample_size / N.
 *
 * @tparam T  Value type to sample.
 */
template <typename T = uint64_t>
struct ReservoirSample {

private:

	std::default_random_engine gen;
    std::vector<T>             samples;
    size_t                     n_sampled;

public:

    ReservoirSample(size_t sample_size)
        : gen       (std::random_device()()),
 		  samples   (sample_size, 0),
          n_sampled (0)
    {
        
    }

    void sample(T item) {
        if (n_sampled < samples.size()) {
            samples[n_sampled] = item;
            ++n_sampled;
            return;
        }
        std::uniform_int_distribution<size_t> U(0, ++n_sampled);
        size_t j = U(gen);
        if (j < samples.size()) {
            samples[j] = item;
        }
    }

    std::vector<T> get_result() const {
        return samples;
    }

    size_t get_n_sampled() const {
        return n_sampled;
    }

    size_t get_sample_size() const {
        return samples.size();
    }

    void clear() {
        std::fill(samples.begin(), samples.end(), 0);
        n_sampled = 0;
    }

};


/**
 * @Synopsis  Named atomic counter.
 */
struct Gauge : public std::atomic_uint64_t {
    const std::string    family;
    const std::string    name;

    Gauge(const std::string& family,
          const std::string& name)
        : family(family),
          name(name)
    {
    }    
};


}
}

#endif
