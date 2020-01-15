/* boink.hh
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "boink/hashing/hashshifter.hh"

    /*
template<class Derived,
         class HashType,
         class Alphabet>
HashType
HashShifter<Derived, HashType, Alphabet>::set_cursor(const std::string& sequence) {
    if (sequence.length() < _K) {
        throw SequenceLengthException("Sequence must at least length K");
    }
    if (!initialized) {
        load(sequence);
        derived().init();
    } else {
        for (auto cit = sequence.begin(); cit < sequence.begin() + _K; ++cit) {
            shift_right(*cit);
        }
    }
    return get();
}

template<class Derived,
         class HashType,
         class Alphabet>
HashType
HashShifter<Derived, HashType, Alphabet>::set_cursor(const char * sequence) {
    // less safe! does not check length
    if(!initialized) {
        load(sequence);
        derived().init();
    } else {
        for (uint16_t i = 0; i < this->_K; ++i) {
            shift_right(sequence[i]);
        }
    }
    return get();
}


template <class Derived,
          class HashType,
          class Alphabet>
template <typename Iterator>
HashType
HashShifter<Derived, HashType, Alphabet>::set_cursor(const Iterator begin, const Iterator end) {
    Iterator _begin = begin, _end = end;
    if(!initialized) {
        load(begin, end);
        derived().init();
    } else {
        uint16_t l = 0;
        while (_begin != _end) {
            shift_right(*_begin);
            ++_begin;
            ++l;
        }
        if (l < this->_K) {
            throw SequenceLengthException("Sequence must at least length K");
        }
    }
    return get();
}
*/

