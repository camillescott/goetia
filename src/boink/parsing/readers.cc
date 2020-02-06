/**
 * (c) Camille Scott, 2019
 * File   : readers.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 21.10.2019
 */

#include "boink/parsing/readers.hh"
#include "boink/sequences/alphabets.hh"

namespace boink::parsing {


template<class Alphabet>
FastxParser<Alphabet>::FastxParser(const std::string& infile, bool strict) 
    : _filename(infile),
        _spin_lock(0),
        _n_parsed(0),
        _have_qualities(false),
        _is_complete(false),
        _strict(strict),
        _n_skipped(0)
{
    _fp = gzopen(_filename.c_str(), "r");
    _kseq = kseq_init(_fp);

    __asm__ __volatile__ ("" ::: "memory");
}

template<class Alphabet>
FastxParser<Alphabet>::~FastxParser() {
    kseq_destroy(_kseq);
    gzclose(_fp);
}

    template class FastxParser<DNA_SIMPLE>;
    template class FastxParser<DNAN_SIMPLE>;
    template class FastxParser<IUPAC_NUCL>;

    template class SplitPairedReader<parsing::FastxParser<DNA_SIMPLE>>;
    template class SplitPairedReader<parsing::FastxParser<DNAN_SIMPLE>>;
    template class SplitPairedReader<parsing::FastxParser<IUPAC_NUCL>>;

}
