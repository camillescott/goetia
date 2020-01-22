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

namespace boink::hashing {

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
struct HashModel {
    
    typedef ValueType value_type;

    value_type hash;

    HashModel(value_type hash)
        : hash(hash)
    {
    }

    HashModel()
        : hash(std::numeric_limits<value_type>::max())
    {
    }

    const value_type value() const {
        return hash;
    }

    operator value_type() const {
        return value();
    }

    __attribute__((visibility("default")))
    friend bool operator>(const HashModel& lhs, const HashModel& rhs) {
        return lhs.value() > rhs.value();
    }

    __attribute__((visibility("default")))
    friend bool operator<(const HashModel& lhs, const HashModel& rhs) {
        return lhs.value() < rhs.value();
    }

    __attribute__((visibility("default")))
    friend bool operator==(const HashModel& lhs, const HashModel& rhs) {
        return lhs.value() == rhs.value();
    }
};

template<class ValueType>
__attribute__((visibility("default")))
inline std::ostream&
operator<<(std::ostream& os, const HashModel<ValueType>& _this) {
    os << "<HashModel h=" << _this.value()
       << ">";
    return os;
}


/**
 * @Synopsis  Model for a canonical k-mer. Stores two value_type
 *            as fw and rc hashes; value() provides the canonical
 *            hash. Matches HashModel interface.
 *
 * @tparam ValueType
 */
template<class ValueType>
struct CanonicalModel {

    typedef ValueType value_type;

    value_type fw_hash, rc_hash;

    CanonicalModel()
        : fw_hash(),
          rc_hash()
    {
    }

    CanonicalModel(value_type fw,
                   value_type rc)
        : fw_hash(fw),
          rc_hash(rc)
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
    friend bool operator>(const CanonicalModel& lhs, const CanonicalModel& rhs) {
        return lhs.value() > rhs.value();
    }

    __attribute__((visibility("default")))
    friend bool operator<(const CanonicalModel& lhs, const CanonicalModel& rhs) {
        return lhs.value() < rhs.value();
    }

    __attribute__((visibility("default")))
    friend bool operator==(const CanonicalModel& lhs, const CanonicalModel& rhs) {
        return lhs.value() == rhs.value();
    }

};

template<class ValueType> __attribute__((visibility("default")))
inline std::ostream&
operator<<(std::ostream& os, const CanonicalModel<ValueType>& _this) {
    os << "<CanonicalModel"
       << " fw=" << _this.fw_hash
       << " rc=" << _this.rc_hash
       << " sign=" << _this.sign()
       << ">";
    return os;
}


/**
 * @Synopsis  Helper struct so we can write, for example, 
 *                WmerModel<HashModel<uint64_t>, Unikmer>
 *            instead of
 *                WmerModel<HashModel, uint64_t, Unikmer>
 *            which makes later composition more general and
 *       safer.
 *
 * @tparam T
 * @tparam MinimizerType
 */
template<class T, class MinimizerType>
struct WmerModel;


/**
 * @Synopsis  A HashModel with an associated minimizer.
 *
 */
template<template<class> class HashType, class ValueType,
         class MinimizerType>
struct WmerModel<HashType<ValueType>, MinimizerType> {
   
    typedef HashType<ValueType> hash_type;
    typedef ValueType           value_type;
    typedef MinimizerType       minimizer_type;

    hash_type     hash;
    minimizer_type minimizer;

    WmerModel(const hash_type&     hash,
              const minimizer_type& minimizer)
        : hash(hash),
          minimizer(minimizer)
    {
    }

    WmerModel() {}

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
    friend bool operator>(const WmerModel& lhs, const WmerModel& rhs) {
        return lhs.hash > rhs.hash;
    }

    __attribute__((visibility("default")))
    friend bool operator<(const WmerModel& lhs, const WmerModel& rhs) {
        return lhs.hash < rhs.hash;
    }

    __attribute__((visibility("default")))
    friend bool operator==(const WmerModel& lhs, const WmerModel& rhs) {
        return lhs.hash == rhs.hash;
    }
};


template<template<class> class HashType, class ValueType,
         class MinimizerType>
__attribute__((visibility("default")))
inline std::ostream&
operator<<(std::ostream& os,
           const WmerModel<HashType<ValueType>, MinimizerType>& wmer) {
    os << "<WmerModel"
       << " hash=" << wmer.hash
       << " minimizer=" << wmer.minimizer
       << ">";
    return os;
}


template<typename T>
struct KmerModel;

/**
 * @Synopsis  Models a Kmer: a Hash and its string repr. The hash_type
 *            must comform to the HashModel template interface.
 *
 *            NOTE: The variadic Extras parameter is need so that KmerModel
 *            can accept a WmerModel as its HashType (or for that matter,
 *            anything that comforms to HashModel but takes extra
 *            template parameters). Naturally, Extras can also just
 *            be empty, in which case it matches regular HashModel<uint64_t>
 *            etc.
 */
template<template<class, class...> class HashType, class ValueType, class...Extras>
struct KmerModel<HashType<ValueType, Extras...>> {

    typedef HashType<ValueType, Extras...> hash_type;
    typedef typename hash_type::value_type value_type;

    hash_type hash;
    std::string kmer;

    KmerModel(const hash_type& hash,
              const std::string& kmer)
        : hash(hash),
          kmer(kmer)
    {
    }

	KmerModel()
        : hash(),
          kmer(0, ' ')
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
    friend bool operator>(const KmerModel& lhs, const KmerModel& rhs) {
        return lhs.hash > rhs.hash;
    }

    __attribute__((visibility("default")))
    friend bool operator<(const KmerModel& lhs, const KmerModel& rhs) {
        return lhs.hash < rhs.hash;
    }

    __attribute__((visibility("default")))
    friend bool operator==(const KmerModel& lhs, const KmerModel& rhs) {
        return lhs.hash == rhs.hash;
    }
};


template<template<class, class...> class HashType, class ValueType,
         class...Extras>
__attribute__((visibility("default")))
inline std::ostream&
operator<<(std::ostream& os, const KmerModel<HashType<ValueType, Extras...>>& kmer) {
    os << "<KmerModel"
       << " hash=" << kmer.hash
       << " kmer=" << kmer.kmer
       << ">";
    return os;
}



/**
 * @Synopsis  helper struct so we can write, for example, 
 *                ShiftModel<HashModel<uint64_t>, DIR_LEFT>
 *            instead of
 *                ShiftModel<HashModel, uint64_t, DIR_LEFT>
 *            which makes later composition more general and
 *            safer. See specialization below.
 *
 * @tparam T
 * @tparam Direction 
 */
template <class T, bool Direction>
struct ShiftModel;

/**
 * @Synopsis  Models a directional shift; conforms to HashModel
 *            and provides the direction and extension symbol.
 *
 */
template <template<class, class...> class HashType, class ValueType, class... Extras,
          bool Direction>
struct ShiftModel<HashType<ValueType, Extras...>, Direction> {
    
    typedef HashType<ValueType, Extras...> hash_type;
    typedef typename hash_type::value_type value_type;
    static const bool direction = Direction;

    hash_type hash;
    char symbol;

    ShiftModel(const hash_type hash,
               const char symbol)
        : hash(hash),
          symbol(symbol)
    {
    }

    ShiftModel()
        : symbol('\0')
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
    friend bool operator>(const ShiftModel& lhs, const ShiftModel& rhs) {
        return lhs.hash > rhs.hash;
    }

    __attribute__((visibility("default")))
    friend bool operator<(const ShiftModel& lhs, const ShiftModel& rhs) {
        return lhs.hash < rhs.hash;
    }

    __attribute__((visibility("default")))
    friend bool operator==(const ShiftModel& lhs, const ShiftModel& rhs) {
        return lhs.hash == rhs.hash;
    }
};

template <template<class, class...> class HashType, class ValueType, class... Extras,
          bool Direction>
__attribute__((visibility("default")))
inline std::ostream&
operator<<(std::ostream& os,
           const ShiftModel<HashType<ValueType, Extras...>, Direction>& shift) {
    os << "<ShiftModel"
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

    Partitioned(value_type& v,
                uint64_t   partition)
        : v(v),
          partition(partition)
    {
    }
    
    Partitioned()
        : partition(std::numeric_limits<uint64_t>::max())
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


typedef boink::hashing::HashModel<uint64_t> Hash;
typedef boink::hashing::CanonicalModel<uint64_t> Canonical;

typedef boink::hashing::KmerModel<boink::hashing::HashModel<uint64_t>> Kmer;
typedef boink::hashing::KmerModel<boink::hashing::CanonicalModel<uint64_t>> CanonicalKmer;

typedef boink::hashing::Partitioned<boink::hashing::Hash> Unikmer;
typedef boink::hashing::Partitioned<boink::hashing::Canonical> CanonicalUnikmer;

typedef boink::hashing::WmerModel<boink::hashing::Hash, boink::hashing::Unikmer> UnikmerWmer;
typedef boink::hashing::WmerModel<boink::hashing::Canonical, boink::hashing::CanonicalUnikmer> CanonicalUnikmerWmer;

typedef boink::hashing::ShiftModel<boink::hashing::Hash, boink::hashing::DIR_LEFT> LeftShift;
typedef boink::hashing::ShiftModel<boink::hashing::Hash, boink::hashing::DIR_RIGHT> RightShift;
typedef boink::hashing::ShiftModel<boink::hashing::Canonical, boink::hashing::DIR_LEFT> LeftCanonicalShift;
typedef boink::hashing::ShiftModel<boink::hashing::Canonical, boink::hashing::DIR_RIGHT> RightCanonicalShift;

}

#endif
