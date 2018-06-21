/* boink.hh
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_HH
#define BOINK_HH

#include <exception>
#include <string>
#include <iostream>
#include <vector>
#include <type_traits>
#include <iterator>
#include <set>

namespace boink {

#define _cerr(x) do { std::ostringstream stream; \
                      stream << x << std::endl; \
                      std::cerr << stream.str(); \
                    } while(0)
#define _cout(x) do { std::ostringstream stream; \
                      stream << x << std::endl; \
                      std::cout << stream.str(); \
                    } while(0)

using std::pair;
using std::vector;


typedef uint64_t hash_t;
typedef std::pair<hash_t, hash_t> full_hash_t;

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

enum direction_t {
    DIR_LEFT,
    DIR_RIGHT
};

typedef vector<hash_t> HashVector;

typedef uint64_t id_t;
typedef pair<hash_t, hash_t> junction_t;


// hash_combine and pair_hash courtesy SO:
// https://stackoverflow.com/questions/9729390/how-to-use-
// unordered-set-with-custom-types/9729747#9729747
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
}


struct junction_hash
{
    inline std::size_t operator()(const junction_t & v) const
    {
         std::size_t seed = 0;
         hash_combine(seed, v.first);
         hash_combine(seed, v.second);
         return seed;
    }
};

std::string junction_to_string(const junction_t& junc) {
    return "(" + std::to_string(junc.first) + ","
               + std::to_string(junc.second) + ")";
}

template<typename _Ty1, typename _Ty2>
std::ostream& operator<<(std::ostream& _os, const std::pair<_Ty1, _Ty2>& _p) {
    _os << "(" << _p.first << ", " << _p.second << ")";
    return _os;
}




class BoinkException : public std::exception {
public:
    explicit BoinkException(const std::string& msg = "Generic boink exception")
        : _msg(msg) { }

    virtual ~BoinkException() throw() { }
    virtual const char* what() const throw ()
    {
        return _msg.c_str();
    }

protected:
    const std::string _msg;
};

} // namespace boink


#endif
