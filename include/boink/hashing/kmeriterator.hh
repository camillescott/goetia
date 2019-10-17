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

#ifndef BOINK_KMERITERATOR_HH
#define BOINK_KMERITERATOR_HH

#include <string>

#include "boink/boink.hh"
#include "boink/hashing/hashshifter.hh"
#include "boink/hashing/exceptions.hh"
#include "boink/kmers/kmerclient.hh"


namespace boink {
namespace hashing {

template <class ShifterType>
class KmerIterator : public kmers::KmerClient {

    const std::string _seq;
    unsigned int index;
    unsigned int length;
    bool _initialized, _shifter_owner;

public:

    typedef ShifterType shifter_type;
    typedef typename ShifterType::hash_type hash_type;

    ShifterType * shifter;

    template<typename... Args>
    KmerIterator(const std::string& seq, uint16_t K, Args&&... args)
        : KmerClient(K), 
          _seq(seq), 
          index(0), 
          _initialized(false), 
          _shifter_owner(true)
    {

        if (_seq.length() < _K) {
            throw SequenceLengthException("Sequence must have length >= K");
        }
        shifter = new ShifterType(seq, K, std::forward<Args>(args)...);
    }

    KmerIterator(const std::string& seq, ShifterType * shifter);
    KmerIterator(const std::string& seq, ShifterType& shifter_proto);

    ~KmerIterator() {
        if (_shifter_owner) {
            delete shifter;
        }
    }

    hash_type first();

    hash_type next();

    bool done() const;

    unsigned int get_start_pos() const {
        if (!_initialized) { return 0; }
        return index - 1;
    }

    unsigned int get_end_pos() const {
        if (!_initialized) { return _K; }
        return index + _K - 1;
    }
};


/*
class FullRollingHasher {
    const char * _seq;
    const std::string _rev;
    const char _ksize;
    unsigned int index;
    unsigned int length;
    bool _initialized;
    oxli::CyclicHash<uint64_t> fwd_hasher;
    oxli::CyclicHash<uint64_t> bwd_hasher;

public:
    FullRollingHasher(const char * seq, unsigned char k) :
        _seq(seq), _rev(oxli::_revcomp(seq)), _ksize(k), index(0),
        _initialized(false), fwd_hasher(k), bwd_hasher(k)
    {
        length = strlen(_seq);
    };

    full_hash_t first() {
        _initialized = true;

        for (char i = 0; i < _ksize; ++i) {
            fwd_hasher.eat(*(_seq + i));
            bwd_hasher.eat(_rev[length - _ksize + i]);
        }
        index += 1;
        return std::make_pair(fwd_hasher.hashvalue, bwd_hasher.hashvalue);
    }

    full_hash_t next() {
        if (!_initialized) {
            return first();
        }

        if (done()) {
            throw oxli_exception("past end of iterator");
        }
        fwd_hasher.update(*(_seq + index - 1), *(_seq + index + _ksize - 1));

        // first argument is added, second is removed from the hash
        bwd_hasher.reverse_update(
          _rev[length - _ksize - index], _rev[length - index]);

        index += 1;

        return std::make_pair(fwd_hasher.hashvalue, bwd_hasher.hashvalue);
    }

    bool done() const {
        return (index + _ksize > length);
    }

    unsigned int get_start_pos() const {
        if (!_initialized) { return 0; }
        return index - 1;
    }

    unsigned int get_end_pos() const {
        if (!_initialized) { return _ksize; }
        return index + _ksize - 1;
    }
};
*/

}
}

#undef pdebug
#endif
