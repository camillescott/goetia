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
#include <iterator>
#include <random>

namespace boink {
namespace metrics {

// From https://gist.github.com/cbsmith/5538174
template <typename RandomGenerator = std::default_random_engine>
struct random_selector
{
    //On most platforms, you probably want to use std::random_device("/dev/urandom")()
    random_selector(RandomGenerator g = RandomGenerator(std::random_device()()))
        : gen(g) {}

    template <typename Iter>
    Iter select(Iter start, Iter end) {
        std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
        std::advance(start, dis(gen));
        return start;
    }

    //convenience function
    template <typename Iter>
    Iter operator()(Iter start, Iter end) {
        return select(start, end);
    }

    //convenience function that works on anything with a sensible begin() and end(), and returns with a ref to the value type
    template <typename Container>
    auto operator()(const Container& c) -> decltype(*begin(c))& {
        return *select(begin(c), end(c));
    }

private:
    RandomGenerator gen;
};


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
    std::vector<T> samples;
    size_t n_sampled;

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


}
}

#endif
