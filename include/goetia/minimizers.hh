/**
 * (c) Camille Scott, 2019
 * File   : minimizers.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 31.08.2019
 */

#ifndef MINIMIZERS_HH
#define MINIMIZERS_HH

#include <cstdint>
#include <deque>
#include <utility>
#include <vector>

#include "goetia/hashing/hashextender.hh"
#include "goetia/hashing/kmeriterator.hh"
#include "goetia/hashing/rollinghashshifter.hh"
#include "goetia/hashing/canonical.hh"

#include "goetia/parsing/parsing.hh"
#include "goetia/parsing/readers.hh"
#include "goetia/processors.hh"


namespace goetia {


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


template <class T = uint64_t>
class InteriorMinimizer : public RollingMin<T> {

protected:

    std::vector<typename RollingMin<T>::value_type> minimizers;

public:

    using typename RollingMin<T>::value_type;
    using RollingMin<T>::RollingMin;
    typedef std::vector<typename RollingMin<T>::value_type> vector_type;

    std::pair<T, int64_t> update(T new_value) {
        auto current = RollingMin<T>::update(new_value);
        if (this->_current_index >= this->_window_size) {
            //if (minimizers.empty() ||
            //    current != minimizers.back()) {

                minimizers.push_back(current);
            //}
        }
        return current;
    }

    const size_t size() const {
        return minimizers.size();
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

    T get_front_value() const {
        return minimizers.front().first;
    }

    T get_back_value() const {
        return minimizers.back().first;
    }

    value_type get_front_minimizer() const {
        return minimizers.front();
    }

    value_type get_back_minimizer() const {
        return minimizers.back();
    }

    void reset() {
        RollingMin<T>::reset();
        minimizers.clear();
    }
};


template <class ShifterType>
struct WKMinimizer {

    typedef ShifterType                       shifter_type;
    _goetia_model_typedefs_from_shiftertype(shifter_type)

    typedef InteriorMinimizer<value_type> minimizer_type;

    class Minimizer : public minimizer_type {

    public:

        const uint16_t K;

        Minimizer(int64_t window_size,
                  uint16_t K)
            : minimizer_type(window_size),
              K(K) {
        }

        auto get_minimizers(const std::string& sequence)
        -> typename minimizer_type::vector_type {
            KmerIterator<ShifterType> iter(sequence, K);
            this->reset();
            
            while(!iter.done()) {
                hash_type h = iter.next();
                minimizer_type::update(h.value());
            }

            return minimizer_type::get_minimizers();
        }

        std::vector<value_type> get_minimizer_values(const std::string& sequence) {
            KmerIterator<ShifterType> iter(sequence, K);
            this->reset();

            while(!iter.done()) {
                hash_type h = iter.next();
                minimizer_type::update(h.value());
            }

            return minimizer_type::get_minimizer_values();
        }

        std::vector<std::pair<std::string, int64_t>> get_minimizer_kmers(const std::string& sequence) {
            typename minimizer_type::vector_type minimizers = get_minimizers(sequence);
            std::vector<std::pair<std::string, int64_t>> kmers;
            for (auto min : minimizers) {
                kmers.push_back(std::make_pair(sequence.substr(min.second, K),
                                               min.second));
            }
            return kmers;
        }
    };


    class Processor : public FileProcessor<Processor,
                                           FastxParser<>> {

    protected:

        Minimizer      M;
        std::string   _output_filename;
        std::ofstream _output_stream;

        typedef FileProcessor<Processor, FastxParser<>> Base;
    public:

        using Base::process_sequence;

        Processor(int32_t window_size,
                  uint16_t K,
                  const std::string& output_filename,
                  uint64_t interval = IntervalCounter::DEFAULT_INTERVAL)
            : Base(interval),
              M(window_size, K),
              _output_filename(output_filename),
              _output_stream(_output_filename.c_str()) {
        
        }

        ~Processor() {
            _output_stream.close();
        }

        uint64_t process_sequence(const Record& read) {
            auto minimizers = M.get_minimizers(read.sequence);

            for (const auto& min : minimizers) {
                _output_stream << this->n_sequences() << ","
                               << min.second << ","
                               << min.first << ","
                               << read.sequence.substr(min.second, M.K)
                               << std::endl;
            }

            return minimizers.size();
        }

        void report() {}
    };

};

extern template class InteriorMinimizer<uint64_t>;
extern template class WKMinimizer<FwdLemireShifter>;
extern template class WKMinimizer<CanLemireShifter>;

}


#endif
