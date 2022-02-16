/**
 * (c) Camille Scott, 2019
 * File   : hashshifter.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.07.2019
 */

#ifndef GOETIA_HASHSHIFTER_HH
#define GOETIA_HASHSHIFTER_HH

#include <algorithm>
#include <cstring>
#include <deque>
#include <iterator>
#include <string>
#include <vector>

#include "goetia/hashing/canonical.hh"
#include "goetia/goetia.hh"

#include "goetia/sequences/alphabets.hh"
#include "goetia/sequences/exceptions.hh"

#include "goetia/meta.hh"

#include "goetia/hashing/rollinghashshifter.hh"


namespace goetia {


class UninitializedShifterException : public GoetiaException {
public:
    explicit UninitializedShifterException(const std::string& msg = "Shifter used without hash_base being called.")
        : GoetiaException(msg)
    {
    }
};


template<class ShifterType>
struct has_minimizer {
    static const bool value = false;
};


/**
 * @Synopsis  Policy client class for shifter adapters. Implementations
 *            must define _shift_left, _shift_right, _get, _hash_base,
 *            and _hash. This class is mostly stateless; it just does basic
 *            error handling and dispatches to the implementation.
 *
 * @tparam Derived   The implementation type.
 * @tparam HashType  Hash return type. Specialize for canonical.
 * @tparam Alphabet  The alphabet to hash over.
 */

template<class T>
struct HashShifter;

template <template<typename, typename> typename ShiftPolicy,
                                       typename HashType,
                                       typename Alphabet>
class HashShifter<ShiftPolicy<HashType, Alphabet>>
    : public ShiftPolicy<HashType, Alphabet>,
      public Tagged<HashShifter<ShiftPolicy<HashType, Alphabet>>> {

protected:

    bool initialized;
    typedef Tagged<HashShifter<ShiftPolicy<HashType, Alphabet>>> tagged_type;

public:

    typedef HashShifter<ShiftPolicy<HashType, Alphabet>> type;
    typedef type                                         shifter_type;

    typedef ShiftPolicy<HashType, Alphabet>              shift_policy;
    typedef typename shift_policy::value_type            value_type;
    typedef typename shift_policy::hash_type             hash_type;
    typedef typename shift_policy::kmer_type             kmer_type;
    typedef Alphabet                                     alphabet;
    static constexpr bool has_kmer_span = shift_policy::has_kmer_span;

    using tagged_type::NAME;
    using tagged_type::OBJECT_ABI_VERSION;

    using shift_policy::K;

    friend shift_policy;

    template<typename... ExtraArgs>
    explicit HashShifter(const std::string& start,
                         uint16_t           K,
                         ExtraArgs&&...     args)
        : shift_policy(K, std::forward<ExtraArgs>(args)...),
          initialized(false)
    {
        hash_base(start);
    }

    template<typename... ExtraArgs>
    explicit HashShifter(uint16_t        K,
                         ExtraArgs&&... args)
        : shift_policy(K, std::forward<ExtraArgs>(args)...),
          initialized(false)
    {
    }

    explicit HashShifter(const HashShifter& other)
        : ShiftPolicy<HashType, Alphabet>(static_cast<const shift_policy &>(other)),
          initialized(false)
    {
    }

    HashShifter() = delete;

    hash_type get() {
        //std::cout << "HashShifter get()" << std::endl;
        return this->get_impl();
    }

    hash_type shift_right(const char& out, const char& in) {
        if (!initialized) {
            throw UninitializedShifterException();
        }
        return this->shift_right_impl(out, in);
    }

    hash_type shift_left(const char& in, const char& out) {
        if (!initialized) {
            throw UninitializedShifterException();
        }
        return this->shift_left_impl(in, out);
    }

    hash_type hash_base(const std::string& sequence) {
        //if (sequence.length() < K) {
        //    throw SequenceLengthException("HashShifter::hash_base: Sequence must at least length K");
        //}
        assert(sequence.length() >= this->K);

        hash_type h = hash_base(sequence.c_str());
        return h;
    }

    template<class It>
    hash_type hash_base(It begin, It end) {
        //if (std::distance(begin, end) != K) {
        //    throw SequenceLengthException("HashShifter::hash_base: Iterator distance must be length K");
        //}
        assert(std::distance(begin, end) == K);
        hash_type h = this->hash_base_impl(begin, end);
        initialized = true;
        return h;
    }

    hash_type hash_base(const char * sequence) {
        auto h = this->hash_base_impl(sequence);
        initialized = true;
        return h;
    }

    const hash_type hash(const std::string& sequence) const {
        type hasher(*this);
        return hasher.hash_base(sequence);
    }

    template<typename... ExtraArgs>
    static hash_type hash(const std::string& sequence,
                          const uint16_t K,
                          ExtraArgs&&... args) {
        if (sequence.length() < K) {
            throw SequenceLengthException("HashShifter::hash: Sequence must at least length K");
        }
        //std::cout << "static HashShifter::hash" << std::endl;
        return hash(sequence.c_str(), K, std::forward<ExtraArgs>(args)...);
    }

    template<typename... ExtraArgs>
    static hash_type hash(const char * sequence,
                          const uint16_t K,
                          ExtraArgs&&... args) {
        
        type hasher(K, std::forward<ExtraArgs>(args)...);
        return hasher.hash_base(sequence);
    }

    bool is_initialized() const {
        return initialized;
    }
};

extern template class HashShifter<FwdLemirePolicy>;
extern template class HashShifter<CanLemirePolicy>;

typedef HashShifter<FwdLemirePolicy> FwdLemireShifter;
typedef HashShifter<CanLemirePolicy> CanLemireShifter;

} // goetia

#endif
