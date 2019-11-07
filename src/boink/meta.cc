/**
 * (c) Camille Scott, 2019
 * File   : meta.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 05.11.2019
 */

#include "boink/meta.hh"

const unsigned int boink::Versioned::VERSION = 0;
template <typename T> const std::string boink::Tagged<T>::NAME = "BoinkClass";

