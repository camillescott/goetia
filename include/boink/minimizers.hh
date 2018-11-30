/* minimizers.hh -- rolling basic minimizers
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef MINIMIZERS_HH
#define MINIMIZERS_HH

#include <deque>
#include <utility>
#include <vector>

#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/kmeriterator.hh"


namespace boink {


template <class T>
class RollingMin {

protected:

    std::deque<T> values;
    std::deque<int64_t> indices;
    const int64_t _window_size;
    int64_t  _current_index;

public:

    typedef std::pair<T, int64_t> value_type;

    RollingMin(int64_t window_size)
        : _window_size(window_size),
          _current_index(0) {

    }

    void reset() {
        values.clear();
        indices.clear();
        _current_index = 0;
    }

    const int64_t window_size() const {
        return _window_size;
    }

    value_type update(T new_value) {

        if (!indices.empty() && 
            (indices.front() <= _current_index - _window_size)) {
            
            values.pop_front();
            indices.pop_front();
        }

        while (!values.empty() &&
               (values.back() > new_value)) {

            values.pop_back();
            indices.pop_back();
        }

        indices.push_back(_current_index);
        values.push_back(new_value);
        ++_current_index;

        return std::make_pair(values.front(), indices.front());
    }

};


template <class T>
class InteriorMinimizer : public RollingMin<T> {

protected:

    std::vector<typename RollingMin<T>::value_type> minimizers;

public:

    using typename RollingMin<T>::value_type;
    using RollingMin<T>::RollingMin;

    std::pair<T, int64_t> update(T new_value) {
        auto current = RollingMin<T>::update(new_value);
        if (this->_current_index >= this->_window_size) {
            if (minimizers.empty() ||
                current != minimizers.back()) {
                minimizers.push_back(current);
            }
        }
        return current;
    }

    std::vector<value_type> get_minimizers() const {
        return minimizers;
    }

    std::vector<T> get_minimizer_values() const {
        std::vector<T> values;
        for (auto m : minimizers) {
            values.push_back(m.first);
        }
        return values;
    }

    void reset() {
        RollingMin<T>::reset();
        minimizers.clear();
    }
};


template <class ShifterType>
class WKMinimizer : public InteriorMinimizer<hashing::hash_t>,
                    public hashing::KmerClient {

public:

    using InteriorMinimizer<hashing::hash_t>::InteriorMinimizer;
    using typename InteriorMinimizer<hashing::hash_t>::value_type;

    WKMinimizer(int64_t window_size,
                uint16_t K)
        : InteriorMinimizer<hashing::hash_t>(window_size),
          hashing::KmerClient(K) {
    }

    std::vector<value_type> get_minimizers(const std::string& sequence) {
        hashing::KmerIterator<ShifterType> iter(sequence, this->_K);
        this->reset();
        
        while(!iter.done()) {
            hashing::hash_t h = iter.next();
            InteriorMinimizer<hashing::hash_t>::update(h);
        }

        return InteriorMinimizer<hashing::hash_t>::get_minimizers();
    }

    std::vector<hashing::hash_t> get_minimizer_values(const std::string& sequence) {
        hashing::KmerIterator<ShifterType> iter(sequence, this->_K);
        this->reset();

        while(!iter.done()) {
            hashing::hash_t h = iter.next();
            InteriorMinimizer<hashing::hash_t>::update(h);
        }

        return InteriorMinimizer<hashing::hash_t>::get_minimizer_values();
    }

    std::vector<std::pair<std::string, int64_t>> get_minimizer_kmers(const std::string& sequence) {
        std::vector<value_type> minimizers = get_minimizers(sequence);
        std::vector<std::pair<std::string, int64_t>> kmers;
        for (auto min : minimizers) {
            kmers.push_back(std::make_pair(sequence.substr(min.second, this->_K),
                                           min.second));
        }
        return kmers;
    }
};


}


#endif
