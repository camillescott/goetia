/* hashing_types.hh -- k-mer hash functions
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_HASHING_TYPES_HH
#define BOINK_HASHING_TYPES_HH

#include <iostream>

namespace boink {
namespace hashing {

# ifdef DEBUG_HASHING
#   define pdebug(x) do { std::cerr << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif

// Type for hashing k-mers into
typedef uint64_t hash_t;
// Type for storing forward and rc hashes of a k-mer
typedef std::pair<hash_t, hash_t> full_hash_t;

// Type for representing a neighbor hash with its prefix or suffix symbol
struct shift_t {
    hash_t hash;
    char symbol;

    shift_t() : 
        hash(0),
        symbol('A') {
    
    }

    shift_t(hash_t hash, char symbol) : 
        hash(hash),
        symbol(symbol) {
    
    }
};


std::ostream& operator<<(std::ostream& os, const shift_t& shift)
{
    os << "<shift_t symbol=" << shift.symbol << " hash=" << shift.hash << ">";
    return os;
}

// A k-mer string and its hash value.
struct kmer_t {
    const hash_t hash;
    const std::string kmer;

    kmer_t(const hash_t hash, const std::string kmer) :
        hash(hash),
        kmer(kmer) {
    }

    bool operator==(const kmer_t& other) const {
        return other.hash == this->hash;
    }

};


std::ostream& operator<<(std::ostream& os, const kmer_t& kmer)
{
    os << "<kmer_t kmer=" << kmer.kmer << " hash=" << kmer.hash << ">";
    return os;
}

} // hashing
} // boink

#endif
