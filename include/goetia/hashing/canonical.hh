/**
 * (c) Camille Scott, 2019
 * File   : canonical.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 08.01.2020
 */

#ifndef BOINK_CANONICAL_HH
#define BOINK_CANONICAL_HH

#include <cstdint>
#include <iostream>
#include <limits>

namespace goetia::hashing {

inline constexpr bool DIR_LEFT = false;
inline constexpr bool DIR_RIGHT = true;


/**
 * @Synopsis  Default model for a regular, single-value hash.
 *            Generally, just a thin wrapper for uint64_t
 *            which provides a common interface exposing
 *            the underlying value through value().
 *
 * @tparam ValueType Underlying storage type.
 */
template<class ValueType = uint64_t>
struct Hash {
    
    typedef ValueType value_type;

    value_type hash;

    Hash(const value_type hash)
        : hash(hash)
    {
    }

    Hash()
        : hash(std::numeric_limits<value_type>::max())
    {
    }

    Hash(const Hash& other)
        : hash(other.hash)
    {
    }

    const value_type value() const {
        return hash;
    }

    operator value_type() const {
        return value();
    }

    __attribute__((visibility("default")))
    friend bool operator>(const Hash& lhs, const Hash& rhs) {
        return lhs.value() > rhs.value();
    }

    __attribute__((visibility("default")))
    friend bool operator<(const Hash& lhs, const Hash& rhs) {
        return lhs.value() < rhs.value();
    }

    __attribute__((visibility("default")))
    friend bool operator==(const Hash& lhs, const Hash& rhs) {
        return lhs.value() == rhs.value();
    }
};

template<class ValueType>
__attribute__((visibility("default")))
inline std::ostream&
operator<<(std::ostream& os, const Hash<ValueType>& _this) {
    os << "<Hash h=" << _this.value()
       << ">";
    return os;
}


/**
 * @Synopsis  Model for a canonical k-mer. Stores two value_type
 *            as fw and rc hashes; value() provides the canonical
 *            hash. Matches Hash interface.
 *
 * @tparam ValueType
 */
template<class ValueType = uint64_t>
struct Canonical {

    typedef ValueType value_type;

    value_type fw_hash, rc_hash;

    Canonical()
        : fw_hash(std::numeric_limits<value_type>::max()),
          rc_hash(std::numeric_limits<value_type>::max())
    {
    }

    Canonical(const value_type fw,
                   const value_type rc)
        : fw_hash(fw),
          rc_hash(rc)
    {
    }

    Canonical(const Canonical& other)
        : fw_hash(other.fw_hash),
          rc_hash(other.rc_hash)
    {
    }

    const bool sign() const {
        return fw_hash < rc_hash ? true : false;
    }

    const value_type value() const {
        return fw_hash < rc_hash ? fw_hash : rc_hash;
    }

    operator value_type() const {
        return value();
    }

    __attribute__((visibility("default")))
    friend bool operator>(const Canonical& lhs, const Canonical& rhs) {
        return lhs.value() > rhs.value();
    }

    __attribute__((visibility("default")))
    friend bool operator<(const Canonical& lhs, const Canonical& rhs) {
        return lhs.value() < rhs.value();
    }

    __attribute__((visibility("default")))
    friend bool operator==(const Canonical& lhs, const Canonical& rhs) {
        return lhs.value() == rhs.value();
    }

};

template<class ValueType> __attribute__((visibility("default")))
inline std::ostream&
operator<<(std::ostream& os, const Canonical<ValueType>& _this) {
    os << "<Canonical"
       << " fw=" << _this.fw_hash
       << " rc=" << _this.rc_hash
       << " sign=" << _this.sign()
       << ">";
    return os;
}


/**
 * @Synopsis  Helper struct so we can write, for example, 
 *                Wmer<Hash<uint64_t>, Unikmer>
 *            instead of
 *                Wmer<Hash, uint64_t, Unikmer>
 *            which makes later composition more general and
 *       safer.
 *
 * @tparam T
 * @tparam MinimizerType
 */
template<class T, class MinimizerType>
struct Wmer;


/**
 * @Synopsis  A Hash with an associated minimizer.
 *
 */
template<template<class> class HashType, class ValueType,
         class MinimizerType>
struct Wmer<HashType<ValueType>, MinimizerType> {
   
    typedef HashType<ValueType> hash_type;
    typedef ValueType           value_type;
    typedef MinimizerType       minimizer_type;

    hash_type     hash;
    minimizer_type minimizer;

    Wmer(const hash_type&     hash,
              const minimizer_type& minimizer)
        : hash(hash),
          minimizer(minimizer)
    {
    }

    Wmer()
        : hash(),
          minimizer()
    {
    }

    Wmer(const Wmer& other)
        : hash(other.hash),
          minimizer(other.minimizer)
    {
    }

    const value_type value() const {
        return hash.value();
    }

    operator value_type() const {
        return hash.value();
    }

    operator hash_type() const {
        return hash;
    }

    __attribute__((visibility("default")))
    friend bool operator>(const Wmer& lhs, const Wmer& rhs) {
        return lhs.hash > rhs.hash;
    }

    __attribute__((visibility("default")))
    friend bool operator<(const Wmer& lhs, const Wmer& rhs) {
        return lhs.hash < rhs.hash;
    }

    __attribute__((visibility("default")))
    friend bool operator==(const Wmer& lhs, const Wmer& rhs) {
        return lhs.hash == rhs.hash;
    }
};


template<template<class> class HashType, class ValueType,
         class MinimizerType>
__attribute__((visibility("default")))
inline std::ostream&
operator<<(std::ostream& os,
           const Wmer<HashType<ValueType>, MinimizerType>& wmer) {
    os << "<Wmer"
       << " hash=" << wmer.hash
       << " minimizer=" << wmer.minimizer
       << ">";
    return os;
}


template<typename T>
struct Kmer;

/**
 * @Synopsis  Models a Kmer: a Hash and its string repr. The hash_type
 *            must comform to the Hash template interface.
 *
 *            NOTE: The variadic Extras parameter is need so that Kmer
 *            can accept a Wmer as its HashType (or for that matter,
 *            anything that comforms to Hash but takes extra
 *            template parameters). Naturally, Extras can also just
 *            be empty, in which case it matches regular Hash<uint64_t>
 *            etc.
 */
template<template<class, class...> class HashType, class ValueType, class...Extras>
struct Kmer<HashType<ValueType, Extras...>> {

    typedef HashType<ValueType, Extras...> hash_type;
    typedef typename hash_type::value_type value_type;

    hash_type hash;
    std::string kmer;

    Kmer(const hash_type& hash,
              const std::string& kmer)
        : hash(hash),
          kmer(kmer)
    {
    }

	Kmer()
        : hash(),
          kmer(0, ' ')
    {
    }

    Kmer(const Kmer& other)
        : hash(other.hash),
          kmer(other.kmer)
    {
    }

    const value_type value() const {
        return hash.value();
    }

    operator value_type() const {
        return hash.value();
    }

    operator hash_type() const {
        return hash;
    }

    __attribute__((visibility("default")))
    friend bool operator>(const Kmer& lhs, const Kmer& rhs) {
        return lhs.hash > rhs.hash;
    }

    __attribute__((visibility("default")))
    friend bool operator<(const Kmer& lhs, const Kmer& rhs) {
        return lhs.hash < rhs.hash;
    }

    __attribute__((visibility("default")))
    friend bool operator==(const Kmer& lhs, const Kmer& rhs) {
        return lhs.hash == rhs.hash;
    }
};


template<template<class, class...> class HashType, class ValueType,
         class...Extras>
__attribute__((visibility("default")))
inline std::ostream&
operator<<(std::ostream& os, const Kmer<HashType<ValueType, Extras...>>& kmer) {
    os << "<Kmer"
       << " hash=" << kmer.hash
       << " kmer=" << kmer.kmer
       << ">";
    return os;
}



/**
 * @Synopsis  helper struct so we can write, for example, 
 *                Shift<Hash<uint64_t>, DIR_LEFT>
 *            instead of
 *                Shift<Hash, uint64_t, DIR_LEFT>
 *            which makes later composition more general and
 *            safer. See specialization below.
 *
 * @tparam T
 * @tparam Direction 
 */
template <class T, bool Direction>
struct Shift;

/**
 * @Synopsis  Models a directional shift; conforms to Hash
 *            and provides the direction and extension symbol.
 *
 */
template <template<class, class...> class HashType, class ValueType, class... Extras,
          bool Direction>
struct Shift<HashType<ValueType, Extras...>, Direction> {
    
    typedef HashType<ValueType, Extras...> hash_type;
    typedef typename hash_type::value_type value_type;
    static const bool direction = Direction;

    hash_type hash;
    char symbol;

    Shift(const hash_type hash,
               const char symbol)
        : hash(hash),
          symbol(symbol)
    {
    }

    Shift()
        : hash(),
          symbol('\0')
    {
    }

    Shift(const Shift& other)
        : Shift(other.hash, other.symbol)
    {
    }

    const value_type value() const {
        return hash.value();
    }

    operator value_type() const {
        return hash.value();
    }

    operator hash_type() const {
        return hash;
    }

    __attribute__((visibility("default")))
    friend bool operator>(const Shift& lhs, const Shift& rhs) {
        return lhs.hash > rhs.hash;
    }

    __attribute__((visibility("default")))
    friend bool operator<(const Shift& lhs, const Shift& rhs) {
        return lhs.hash < rhs.hash;
    }

    __attribute__((visibility("default")))
    friend bool operator==(const Shift& lhs, const Shift& rhs) {
        return lhs.hash == rhs.hash;
    }
};

template <template<class, class...> class HashType, class ValueType, class... Extras,
          bool Direction>
__attribute__((visibility("default")))
inline std::ostream&
operator<<(std::ostream& os,
           const Shift<HashType<ValueType, Extras...>, Direction>& shift) {
    os << "<Shift"
       << " hash=" << shift.hash
       << " symbol=" << shift.symbol
       << " direction=" << shift.direction
       << ">";
    return os;
}


/**
 * @Synopsis  Models any value with an associated partition ID.
 *
 * @tparam ValueType
 */
template <class ValueType>
struct Partitioned {

    typedef ValueType value_type;

    value_type v;
    uint64_t   partition;

    Partitioned(const value_type& v,
                const uint64_t   partition)
        : v(v),
          partition(partition)
    {
    }
    
    Partitioned()
        : v(),
          partition(std::numeric_limits<uint64_t>::max())
    {
    }

    Partitioned(const Partitioned& other)
        : v(other.v),
          partition(other.partition)
    {
    }

    const value_type value() const {
        return v;
    }

    __attribute__((visibility("default")))
    friend bool operator>(const Partitioned& lhs, const Partitioned& rhs) {
        return lhs.v > rhs.v;
    }

    __attribute__((visibility("default")))
    friend bool operator<(const Partitioned& lhs, const Partitioned& rhs) {
        return lhs.v < rhs.v;
    }

    __attribute__((visibility("default")))
    friend bool operator==(const Partitioned& lhs, const Partitioned& rhs) {
        return lhs.v == rhs.v;
    }
};


template <class ValueType>
__attribute__((visibility("default")))
inline std::ostream&
operator<<(std::ostream& os, const Partitioned<ValueType>& p) {
    os << "<Partitioned"
       << " value=" << p.v
       << " partition=" << p.partition
       << ">";
    return os;
}



typedef goetia::hashing::Partitioned<goetia::hashing::Hash<>> Unikmer;
typedef goetia::hashing::Partitioned<goetia::hashing::Canonical<>> CanonicalUnikmer;

typedef goetia::hashing::Wmer<goetia::hashing::Hash<>, goetia::hashing::Unikmer> UnikmerWmer;
typedef goetia::hashing::Wmer<goetia::hashing::Canonical<>, goetia::hashing::CanonicalUnikmer> CanonicalUnikmerWmer;

typedef goetia::hashing::Shift<goetia::hashing::Hash<>, goetia::hashing::DIR_LEFT> LeftShift;
typedef goetia::hashing::Shift<goetia::hashing::Hash<>, goetia::hashing::DIR_RIGHT> RightShift;
typedef goetia::hashing::Shift<goetia::hashing::Canonical<>, goetia::hashing::DIR_LEFT> LeftCanonicalShift;
typedef goetia::hashing::Shift<goetia::hashing::Canonical<>, goetia::hashing::DIR_RIGHT> RightCanonicalShift;

}

extern template class goetia::hashing::Hash<uint64_t>;
extern template class goetia::hashing::Canonical<goetia::hashing::Hash<uint64_t>>;

extern template class goetia::hashing::Kmer<goetia::hashing::Hash<uint64_t>>;
extern template class goetia::hashing::Kmer<goetia::hashing::Canonical<uint64_t>>;

extern template class goetia::hashing::Kmer<goetia::hashing::UnikmerWmer>;
extern template class goetia::hashing::Kmer<goetia::hashing::CanonicalUnikmerWmer>;

extern template class goetia::hashing::Wmer<goetia::hashing::Hash<uint64_t>, goetia::hashing::Unikmer>;
extern template class goetia::hashing::Wmer<goetia::hashing::Canonical<uint64_t>, goetia::hashing::CanonicalUnikmer>;

extern template class goetia::hashing::Shift<goetia::hashing::Hash<uint64_t>, goetia::hashing::DIR_LEFT>;
extern template class goetia::hashing::Shift<goetia::hashing::Hash<uint64_t>, goetia::hashing::DIR_RIGHT>;
extern template class goetia::hashing::Shift<goetia::hashing::Canonical<uint64_t>, goetia::hashing::DIR_LEFT>;
extern template class goetia::hashing::Shift<goetia::hashing::Canonical<uint64_t>, goetia::hashing::DIR_RIGHT>;



#endif
