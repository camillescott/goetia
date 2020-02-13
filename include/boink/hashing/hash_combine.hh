/**
 * (c) Camille Scott, 2019
 * File   : hash_combine.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 29.10.2019
 *
 * Generic tuple hashing shamelessly riffed from:
 * https://stackoverflow.com/questions/7110301/generic-hash-for-tuples-in-unordered-map-unordered-set
 *
 */

#ifndef HASH_COMBINE_HH
#define HASH_COMBINE_HH

#include <cstdint>
#include <functional>
#include <utility>

namespace hash_tuple {

    using std::size_t;

    template <typename TT>
    struct hash {
        size_t operator()(TT const& tt) const {
            return std::hash<TT>()(tt);
        }
    };

    template <class T>
    inline void hash_combine(std::size_t& seed, T const& v) {
        seed ^= hash_tuple::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    namespace {
        // Recursive template code derived from Matthieu M.
        template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
        struct HashValueImpl {
            static void apply(size_t& seed, Tuple const& tuple) {
                HashValueImpl<Tuple, Index - 1>::apply(seed, tuple);
                hash_combine(seed, std::get<Index>(tuple));
            }
        };

        template <class Tuple>
        struct HashValueImpl<Tuple, 0> {
            static void apply(size_t& seed, Tuple const& tuple) {
                hash_combine(seed, std::get<0>(tuple));
            }
        };
    }  // namespace

    template <typename... TT>
    struct hash<std::tuple<TT...>> {
        size_t operator()(std::tuple<TT...> const& tt) const {
            size_t seed = 0;
            HashValueImpl<std::tuple<TT...>>::apply(seed, tt);
            return seed;
        }
    };
}  // namespace hash_tuple

#endif
