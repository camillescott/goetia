/**
 * (c) Camille Scott, 2019
 * File   : kmeriterator.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 23.07.2019
 */
/* kmeriterator.hh -- sequence k-mer iterators
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef GOETIA_KMERITERATOR_HH
#define GOETIA_KMERITERATOR_HH

#include <iterator>
#include <string>
#include <string_view>

#include "goetia/goetia.hh"
#include "goetia/hashing/hashshifter.hh"
#include "goetia/sequences/exceptions.hh"

#include "goetia/hashing/shifter_types.hh"


namespace goetia {

template <class ShifterType>
class KmerIterator {

    const std::string _seq;
    unsigned int index;
    unsigned int length;
    bool _initialized, _shifter_owner;

public:

    typedef ShifterType                     shifter_type;
    typedef typename ShifterType::hash_type hash_type;
    const uint16_t                          K;

    ShifterType * shifter;

    template<typename... Args>
    __attribute__((visibility("default")))
    explicit KmerIterator(const std::string& seq, uint16_t K, Args&&... args)
        : K(K), 
          _seq(seq), 
          index(0), 
          _initialized(false), 
          _shifter_owner(true)
    {

        if (_seq.length() < K) {
            throw SequenceLengthException("Sequence must have length >= K");
        }
        shifter = new ShifterType(seq, K, std::forward<Args>(args)...);
    }

    __attribute__((visibility("default")))
    explicit KmerIterator(const std::string& seq, ShifterType * shifter)
        : K(shifter->K), 
          _seq(seq),
          index(0), 
          _initialized(false),
          _shifter_owner(false), 
          shifter(shifter) 
    {
        if (_seq.length() < K) {
            throw SequenceLengthException("Sequence must have length >= K");
        }
        shifter->hash_base(_seq);
    }

    __attribute__((visibility("default")))
    explicit KmerIterator(const std::string& seq, ShifterType& shifter_proto)
        : K(shifter_proto.K), 
          _seq(seq), 
          index(0), 
          _initialized(false), 
          _shifter_owner(true)
    {

        if (_seq.length() < K) {
            throw SequenceLengthException("Sequence must have length >= K");
        }
        shifter = new ShifterType(shifter_proto);
        shifter->hash_base(_seq);
    }

    __attribute__((visibility("default")))
    ~KmerIterator() {
        if (_shifter_owner) {
            delete shifter;
        }
    }

    __attribute__((visibility("default")))
    hash_type first()  {
        _initialized = true;
        index += 1;
        return shifter->get();
    }

    __attribute__((visibility("default")))
    hash_type next() {
        if (!_initialized) {
            return first();
        }

        if (done()) {
            //throw InvalidCharacterException("past end of iterator");
            return hash_type();
        }

        auto ret = shifter->shift_right(_seq[index - 1], _seq[index + K - 1]);
        index += 1;

        return ret;
    }
    
    hash_type get() {
        if (!_initialized) {
            return first();
        }

        if (done()) {
            return hash_type();
        }

        return shifter->get();
    }

    __attribute__((visibility("default")))
    bool done() const  {
        return (index + K > _seq.length());
    }

    __attribute__((visibility("default")))
    unsigned int get_start_pos() const {
        if (!_initialized) { return 0; }
        return index - 1;
    }

    __attribute__((visibility("default")))
    unsigned int get_end_pos() const {
        if (!_initialized) { return K; }
        return index + K - 1;
    }
};


template <typename ShifterType>
struct kmer_iterator {
    using value_type         = typename ShifterType::hash_type;
    using iterator_category = std::input_iterator_tag;
    using difference_type   = int64_t;
    using reference         = const value_type;

    kmer_iterator(ShifterType * shifter, size_t start, const std::string_view sequence)
        : shifter(shifter), sequence(sequence), index(start)
    {
        if (!_out_of_bounds()) {
            shifter->hash_base(sequence.begin() + index, 
                               sequence.begin() + index + shifter->K);
        }
    }

    reference operator*() const {
        if (_out_of_bounds()) {
            return value_type();
        }
        return shifter->get();
    }

    kmer_iterator& operator++() {
        index += 1;
        if (!_out_of_bounds()) {
            shifter->shift_right(sequence[index - 1],
                                 sequence[index + shifter->K - 1]);
        }
        return *this;
    }

    kmer_iterator operator++(int) {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    friend bool operator==(const kmer_iterator& a, const kmer_iterator& b) {
        return (*a == *b) && (a.shifter == b.shifter);
    }

    friend bool operator!=(const kmer_iterator& a, const kmer_iterator& b) {
        return (*a != *b) || (a.shifter != b.shifter);
    }

private:
    ShifterType * shifter;
    const std::string_view sequence;
    size_t index;

    const bool _out_of_bounds() const {
        return index + shifter->K > sequence.length();
    }
};


template <typename ShifterType,
          typename... Args>
struct kmer_iterator_wrapper {
    const std::string_view sequence;
    ShifterType shifter;

    kmer_iterator_wrapper(const std::string_view sequence,
                          uint16_t K,
                          Args&&... args)
        : sequence(sequence),
          shifter(K, std::forward<Args>(args)...) {}

    auto begin() {
        return kmer_iterator(&shifter, 0, sequence);
    }

    auto end() {
        return kmer_iterator(&shifter, sequence.length() - shifter.K + 1, sequence);
    }
};


template <typename ShifterType,
          typename... Args>
constexpr auto hash_sequence(const std::string_view sequence,
                             uint16_t K,
                             Args&&... args) {

    return kmer_iterator_wrapper<ShifterType>{ sequence, K, std::forward<Args>(args)... };
}


void test_kmer_iterator();


extern template class KmerIterator<FwdLemireShifter>;
extern template class KmerIterator<CanLemireShifter>;

extern template class KmerIterator<FwdUnikmerShifter>;
extern template class KmerIterator<CanUnikmerShifter>;


}

#undef pdebug
#endif
