/**
 * (c) Camille Scott, 2019
 * File   : kmeriterator.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 29.08.2019
 */

#include "goetia/hashing/kmeriterator.hh"
#include "goetia/hashing/shifter_types.hh"

namespace goetia {

    template class KmerIterator<FwdLemireShifter>;
    template class KmerIterator<CanLemireShifter>;

    template class KmerIterator<FwdUnikmerShifter>;
    template class KmerIterator<CanUnikmerShifter>;


    void test_kmer_iterator() {
        std::string seq = "AAACCATTCTATATCTCTCTTAAAAACTTCTAATA";
        uint16_t K = 21;

        std::cout << "kmer_iterator" << std::endl;
        for (const auto& h : hash_sequence<FwdLemireShifter>(seq, 21)) {
            std::cout << h << " ";
        } std::cout << std::endl;

        std::cout << "legacy KmerIterator" << std::endl;
        KmerIterator<FwdLemireShifter> it(seq, K);
        while (!it.done()) {
            auto h = it.next();
            std::cout << h << " ";
        } std::cout << std::endl;
    }

}
