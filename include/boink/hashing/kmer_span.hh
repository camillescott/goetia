/**
 * (c) Camille Scott, 2019
 * File   : kmer_span.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 10.01.2020
 */

#ifndef BOINK_KMER_SPAN_HH
#define BOINK_KMER_SPAN_HH

#include <iostream>
#include <type_traits>
#include <deque>
#include <string>

#include "boink/ring_span.hpp"

namespace boink::hashing {

void test();


class KmerSpan {

protected:

    char *                  kmer_buffer;
    bool                    loaded;

public:

    nonstd::ring_span<char> ring;
    const uint16_t          K;

    explicit KmerSpan(uint16_t K)
        : K(K),
          kmer_buffer(new char[K]),
          ring(kmer_buffer, kmer_buffer + K, kmer_buffer, K),
          loaded(false)
    {
    }

    ~KmerSpan() {
        delete [] kmer_buffer;
    }

    const bool is_loaded() const {
        return loaded;
    }

    const char& front() const {
        return ring.front();
    }

    const char& back() const {
        return ring.back();
    }

    __attribute__((visibility("default")))
    void load(const std::string& sequence) {
        for (uint16_t i = 0; i < K; ++i)  {
            ring.push_back(sequence[i]);
        }
        loaded = true;
    }

    __attribute__((visibility("default")))
    void load(const char * sequence) {
        for (uint16_t i = 0; i < K; ++i)  {
            ring.push_back(sequence[i]);
        }
        loaded = true;
    }

    template <typename Iterator>
    __attribute__((visibility("default")))
    void load(Iterator begin, Iterator end) {
        while (begin != end) {
            ring.push_back(*begin);
            ++begin;
        }
    }

    template <typename Iterator>
    __attribute__((visibility("default")))
    void load(Iterator begin) {
        for (uint16_t i = 0; i < K; ++i)  {
            ring.push_back(*begin);
            ++begin;
        }
        loaded = true;
    }

    const std::string to_string() const {
        return std::string(ring.begin(),
                           ring.end());
    }

    void to_deque(std::deque<char>& d) const {

        for (auto symbol : ring) {
            d.push_back(symbol);
        }
    }
};

}


#endif
