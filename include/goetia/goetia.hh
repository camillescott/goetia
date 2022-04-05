/**
 * (c) Camille Scott, 2019
 * File   : goetia.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 13.01.2020
 */


#ifndef GOETIA_MISC_HH
#define GOETIA_MISC_HH

#include <cstdint>
#include <string>
#include <iostream>
#include <iterator>
#include <iostream>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace goetia {


#define _cerr(x) do { std::ostringstream stream; \
                      stream << x << std::endl; \
                      std::cerr << stream.str(); \
                    } while(0)
#define _cout(x) do { std::ostringstream stream; \
                      stream << x << std::endl; \
                      std::cout << stream.str(); \
                    } while(0)


template <typename T>
std::string repr(const T& item) {
    std::ostringstream _os;
    static_cast<std::ostream&>(_os) << item;
    return _os.str();
}


template<typename _Ty1, typename _Ty2>
std::string repr(const std::pair<_Ty1, _Ty2>& _p) {
    std::ostringstream _os;
    _os << "(" << _p.first << ", " << _p.second << ")";
    return _os.str();
}


template <template<typename, typename> class Container, class T, class A>
std::string repr(const Container<T,A>& v) {
    std::ostringstream os;
    os << "[";
    /*
    for (size_t i = 0; i < v.size(); ++i) {
        os << repr(v[i]);
        if (i != v.size() - 1)
            os << ", ";
    }
    */
    size_t i = 0;
    for (const auto& elem : v) {
        ++i;
        os << repr(elem);
        if (i != v.size())
            os << ", ";
    }
    os << "]";
    return os.str();
}


template <template<typename, typename, typename> class Container, class T, class C, class A>
std::string repr(const Container<T,C,A>& v) {
    std::ostringstream os;
    os << "[";
    size_t i = 0;
    for (const auto& elem : v) {
        ++i;
        os << repr(elem);
        if (i != v.size())
            os << ", ";
    }
    os << "]";
    return os.str();
}


template<template<typename, typename> class Container, class T, class A>
bool contains(const Container<T,A>& collection,
              T&                    item) {
    return collection.count(item);
}


template<typename T>
bool contains(std::vector<T> collection,
              T item) {
    for (auto& i : collection) {
        if (i == item) {
            return true;
        }
    }
    return false;
}

/*
 * from https://www.reedbeta.com/blog/python-like-enumerate-in-cpp17/
 */
template <typename T,
          typename TIter = decltype(std::begin(std::declval<T>())),
          typename = decltype(std::end(std::declval<T>()))>
constexpr auto enumerate(T && iterable) {
    struct iterator
    {
        size_t i;
        TIter iter;
        bool operator != (const iterator & other) const { return iter != other.iter; }
        void operator ++ () { ++i; ++iter; }
        auto operator * () const { return std::tie(i, *iter); }
    };

    struct iterable_wrapper
    {
        T iterable;
        auto begin() { return iterator{ 0, std::begin(iterable) }; }
        auto end() { return iterator{ 0, std::end(iterable) }; }
    };

    return iterable_wrapper{ std::forward<T>(iterable) };
}


/*
 * Adapted from make_from_tuple: https://en.cppreference.com/w/cpp/utility/make_from_tuple
 */
namespace detail {
template <class T, class Tuple, std::size_t... I>
std::shared_ptr<T> make_shared_from_tuple_impl( Tuple&& t, std::index_sequence<I...> )
{
  return std::make_shared<T>(std::get<I>(std::forward<Tuple>(t))...);
}
} // namespace detail
 

template <class T, class Tuple>
std::shared_ptr<T> make_shared_from_tuple( Tuple&& t )
{
    return detail::make_shared_from_tuple_impl<T>(std::forward<Tuple>(t),
        std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<Tuple>>>{});
}


template <typename... Args>
auto make_tuple_from_tuple(std::tuple<Args...> model_tuple, Args&&... args) {
    return std::make_tuple(std::forward<Args>(args)...);
}


template <typename T>
struct fail : std::false_type 
{
};

} // goetia


#endif
