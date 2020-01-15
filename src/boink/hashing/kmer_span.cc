/**
 * (c) Camille Scott, 2019
 * File   : kmer_span.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 10.01.2020
 */


#include "boink/hashing/kmer_span.hh"

namespace boink::hashing {

struct HasSpan : KmerSpanMixin<>::type {

    HasSpan(int K)
        : KmerSpanMixin<>::type(K) {}
};

struct NoSpan {

};

struct DoubleSpan : HasSpan,
                    KmerSpanMixin<HasSpan>::type
{
    DoubleSpan(int K)
        : HasSpan(K),
          KmerSpanMixin<HasSpan>::type(K) {}
};



void test() {
    HasSpan s(11);
    DoubleSpan d(11);
    //KmerSpanMixin<HasSpan>::type mx(11);
    //mx.K();
}

}
