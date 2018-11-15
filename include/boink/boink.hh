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


std::ostream& operator<<(std::ostream& os, const shift_t& shift)
{
    os << "<shift_t symbol=" << shift.symbol << " hash=" << shift.hash << ">";
    return os;
}


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



enum direction_t {
    DIR_LEFT,
    DIR_RIGHT
};


typedef vector<hash_t> HashVector;
typedef uint64_t id_t;
typedef pair<hash_t, hash_t> junction_t;


enum node_meta_t {
    FULL,
    TIP,
    ISLAND,
    CIRCULAR,
    LOOP,
    TRIVIAL,
    DECISION
};

inline const char * node_meta_repr(node_meta_t meta) {
    switch(meta) {
        case FULL:
            return "FULL";
        case TIP:
            return "TIP";
        case ISLAND:
            return "ISLAND";
        case CIRCULAR:
            return "CIRCULAR";
        case TRIVIAL:
            return "TRIVIAL";
        case LOOP:
            return "LOOP";
        case DECISION:
            return "DECISION";
        default:
            return "UNKNOWN";
    }
}


enum update_meta_t {
    BUILD_UNODE,
    BUILD_DNODE,
    DELETE_UNODE,
    SPLIT_UNODE,
    EXTEND_UNODE,
    CLIP_UNODE,
    MERGE_UNODES
};


template<typename _Ty1, typename _Ty2>
std::ostream& operator<<(std::ostream& _os, const std::pair<_Ty1, _Ty2>& _p) {
    _os << "(" << _p.first << ", " << _p.second << ")";
    return _os;
}

template <template<typename, typename> class Container, class T, class A>
std::ostream& operator<<(std::ostream& os, const Container<T,A>& v)
{
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i != v.size() - 1)
            os << ", ";
    }
    os << "]\n";
    return os;
}

class BoinkException : public std::exception {
public:
    explicit BoinkException(const std::string& msg = "Generic boink exception.")
        : _msg(msg) { }

    virtual ~BoinkException() throw() { }
    virtual const char* what() const throw ()
    {
        return _msg.c_str();
    }

protected:
    const std::string _msg;
};

class InvalidCharacterException : public std::exception {
public:
    explicit InvalidCharacterException(const std::string& msg = "Invalid character encountered.")
        : _msg(msg) { }

    virtual ~InvalidCharacterException() throw() { }
    virtual const char* what() const throw ()
    {
        return _msg.c_str();
    }

protected:
    const std::string _msg;
};

class SequenceLengthException: public std::exception {
public:
    explicit SequenceLengthException(const std::string& msg = "Sequence was too short.")
        : _msg(msg) { }

    virtual ~SequenceLengthException() throw() { }
    virtual const char* what() const throw ()
    {
        return _msg.c_str();
    }

protected:
    const std::string _msg;
};

} // namespace boink


#endif
