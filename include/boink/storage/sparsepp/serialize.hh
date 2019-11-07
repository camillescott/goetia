/**
 * (c) Camille Scott, 2019
 * File   : serialize.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 05.11.2019
 */

#ifndef BOINK_SPP_SERIALIZE_HH
#define BOINK_SPP_SERIALIZE_HH

#include <fstream>

#include "boink/storage/sparsepp/spp.h"

namespace boink {
namespace storage {

class BaseSppSerializer {
public:
    // serialize basic types to FILE
    // -----------------------------
    template <class T>
    bool operator()(std::ofstream *out, const T& value)
    {
        out->write(reinterpret_cast<const char*>(&value), sizeof(value));
        return true;
    }

    template <class T>
    bool operator()(std::ifstream *in, T* value)
    {
        in->read(reinterpret_cast<char*>(value), sizeof(*value));
        return true;
    }


    // serialize spp::sparse_hash_map<K, V, H, E, A> to FILE
    // -----------------------------------------------------
    template <class K, class V, 
              class H = spp::spp_hash<K>, 
              class E = std::equal_to<K>, 
              class A = SPP_DEFAULT_ALLOCATOR<std::pair<const K, V> > >
    bool operator()(std::ofstream* out, const spp::sparse_hash_map<K, V, H, E, A>& value)
    {
        return const_cast<spp::sparse_hash_map<K, V, H, E, A> &>(value).serialize(*this, out);
    }

    template <class K, class V, 
              class H = spp::spp_hash<K>, 
              class E = std::equal_to<K>, 
              class A = SPP_DEFAULT_ALLOCATOR<std::pair<const K, V> > >
    bool operator()(std::ifstream* in, spp::sparse_hash_map<K, V, H, E, A> *value)
    {
        new (value) spp::sparse_hash_map<K, V, H, E, A>();
        return value->unserialize(*this, in);
    }

    template <class V,
              class H = spp::spp_hash<V>,
              class E = std::equal_to<V>,
              class A = SPP_DEFAULT_ALLOCATOR<V> >
    bool operator()(std::ofstream* out, const spp::sparse_hash_set<V, H, E, A>& value)
    {
        return const_cast<spp::sparse_hash_set<V, H, E, A> &>(value).serialize(*this, out);
    }

    template <class V,
              class H = spp::spp_hash<V>,
              class E = std::equal_to<V>,
              class A = SPP_DEFAULT_ALLOCATOR<V> >
    bool operator()(std::ifstream* in, spp::sparse_hash_set<V, H, E, A>* value)
    {
        new (value) spp::sparse_hash_set<V, H, E, A>();
        return value->unserialize(*this, in);
    }

    // serialize std::pair<const A, B> to FILE - needed for maps
    // ---------------------------------------------------------
    template <class A, class B>
    bool operator()(std::ofstream* out, const std::pair<const A, B>& value)
    {
        return (*this)(out, value.first) && (*this)(out, value.second);
    }

    template <class A, class B>
    bool operator()(std::ifstream* in, std::pair<const A, B> *value)
    {
        return (*this)(in, (A *)&value->first) && (*this)(in, &value->second);
    }
};



}
}

#endif
