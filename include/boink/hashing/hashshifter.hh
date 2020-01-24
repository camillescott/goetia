/**
 * (c) Camille Scott, 2019
 * File   : hashshifter.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.07.2019
 */

#ifndef BOINK_HASHSHIFTER_HH
#define BOINK_HASHSHIFTER_HH

#include <algorithm>
#include <cstring>
#include <deque>
#include <iterator>
#include <string>
#include <vector>

#include "boink/kmers/kmerclient.hh"
#include "boink/hashing/canonical.hh"
#include "boink/boink.hh"

#include "boink/sequences/alphabets.hh"
#include "boink/sequences/exceptions.hh"

#include "boink/ring_span.hpp"
#include "boink/meta.hh"

namespace boink::hashing {


class UninitializedShifterException : public BoinkException {
public:
    explicit UninitializedShifterException(const std::string& msg = "Shifter used without hash_base being called.")
        : BoinkException(msg)
    {
    }
};


template<class ShifterType>
struct has_minimizer {
    static const bool value = false;
};


/**
 * @Synopsis  CRTP base class for shifter adapters. Implementations
 *            must define _shift_left, _shift_right, _get, _hash_base,
 *            and _hash. This class is mostly stateless; it just does basic
 *            error handling and dispatches to the implementation.
 *
 * @tparam Derived   The implementation type.
 * @tparam HashType  Hash return type. Specialize for canonical.
 * @tparam Alphabet  The alphabet to hash over.
 */
template <class Derived,
          class HashType = HashModel<uint64_t>,
          class Alphabet = DNA_SIMPLE>
class HashShifter : public kmers::KmerClient,
                    public Tagged<HashShifter<Derived,
                                              HashType,
                                              Alphabet>> {

protected:

    bool initialized;
    typedef Tagged<HashShifter<Derived, HashType, Alphabet>> tagged_type;

public:

    typedef HashShifter<Derived, HashType, Alphabet> type;
    typedef typename HashType::value_type            value_type;
    typedef HashType                                 hash_type;
    typedef KmerModel<hash_type>                     kmer_type;

    typedef Alphabet alphabet;

    using tagged_type::NAME;
    using tagged_type::OBJECT_ABI_VERSION;

    hash_type get() {
        std::cout << "HashShifter get() of " << this << std::endl;
        return derived()._get();
    }

    hash_type shift_right(const char& out, const char& in) {
        if (!initialized) {
            throw UninitializedShifterException();
        }
        return derived()._shift_right(out, in);
    }

    hash_type shift_left(const char& in, const char& out) {
        if (!initialized) {
            throw UninitializedShifterException();
        }
        return derived()._shift_left(in, out);
    }

    hash_type hash_base(const std::string& sequence) {
        if (sequence.length() < _K) {
            throw SequenceLengthException("Sequence must at least length K");
        }

        hash_type h = hash_base(sequence.c_str());
        return h;
    }

    template<class It>
    hash_type hash_base(It begin, It end) {
        if (std::distance(begin, end) != _K) {
            throw SequenceLengthException("Iterator distance must be length K");
        }
        hash_type h = derived()._hash_base(begin, end);
        initialized = true;
        return h;
    }

    hash_type hash_base(const char * sequence) {
        std::cout << "HashShifter::hash_base(const char *) " << sequence << std::endl;
        auto h = derived()._hash_base(sequence);
        initialized = true;
        return h;
    }

    const hash_type hash(const std::string& sequence) const {
        Derived hasher(derived());
        return hasher.hash_base(sequence);
    }

    template<typename... ExtraArgs>
    static hash_type hash(const std::string& sequence,
                          const uint16_t K,
                          ExtraArgs&&... args) {
        if (sequence.length() < K) {
            throw SequenceLengthException("Sequence must at least length K");
        }
        //std::cout << "static HashShifter::hash" << std::endl;
        return hash(sequence.c_str(), K, std::forward<ExtraArgs>(args)...);
    }

    template<typename... ExtraArgs>
    static hash_type hash(const char * sequence,
                          const uint16_t K,
                          ExtraArgs&&... args) {

        return Derived::_hash(sequence, K, std::forward<ExtraArgs>(args)...);
    }

    bool is_initialized() const {
        return initialized;
    }

private:

    explicit HashShifter(const std::string& start,
                         uint16_t           K)
        : KmerClient(K),
          initialized(false)
    {
        hash_base(start);
        std::cout << "HashShifter ctor " << this << std::endl;
    }

    explicit HashShifter(uint16_t K)
        : KmerClient(K),
          initialized(false)
    {
        std::cout << "HashShifter ctor " << this << std::endl;
    }

    ~HashShifter()
    {
        std::cout << "HashShifter dstor" << this << std::endl;
    }

    friend Derived;

    Derived& derived() {
        return *static_cast<Derived*>(this);
    }

    const Derived& derived() const {
        return *static_cast<const Derived*>(this);
    }
};

} // boink

#endif
