/* metrics.hh -- utilities for metrics gathering
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef GOETIA_UTILITY_METRICS_HH
#define GOETIA_UTILITY_METRICS_HH

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iterator>
#include <random>

namespace goetia {

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
struct Gauge : public std::atomic_int64_t {
    const std::string    family;
    const std::string    name;

    Gauge(const std::string& family,
          const std::string& name)
        : std::atomic_int64_t(0),
          family(family),
          name(name)
    {
    }    
};


/**
 * @Synopsis  Bounded counter. Returns True when interval is
 *            reached and resets to zero.
 */
class IntervalCounter {

private:

    uint64_t _counter;
    uint64_t _total;

public:

    const uint64_t interval;

    static constexpr uint64_t DEFAULT_INTERVAL = 500000;

    IntervalCounter(uint64_t interval = DEFAULT_INTERVAL)
        : interval(interval),
          _counter(0),
          _total(0)
    {
    }

    const uint64_t counter() const {
        return _counter;
    }

    const uint64_t total() const {
        return _total;
    }

    bool poll(uint64_t incr=1) {
        _counter += incr;
        _total += incr;
        if (_counter >= interval) {
            _counter = 0;
            return true;
        } else {
            return false;
        }
    }

    friend inline std::ostream& operator<< (std::ostream& os, const IntervalCounter& counter) {
       os  << "IntervalCounter<interval=" << counter.interval 
           << ", counter=" << counter.counter() << ">";
       return os;
    }
};

extern template class ReservoirSample<uint64_t>;
extern template class ReservoirSample<double>;
extern template class ReservoirSample<float>;
extern template class ReservoirSample<uint32_t>;

}

#endif
