/**
 * (c) Camille Scott, 2019
 * File   : kmer_span.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 10.01.2020
 */


#include "boink/hashing/kmer_span.hh"

namespace boink::hashing {

template
void KmerSpanMixinImpl<true>::load<std::string::iterator>(std::string::iterator begin,
                                                          std::string::iterator end);

template
void KmerSpanMixinImpl<true>::load<std::deque<char>::iterator>(std::deque<char>::iterator begin,
                                                               std::deque<char>::iterator end);

}
