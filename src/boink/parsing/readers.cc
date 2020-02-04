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

    template class FastxParser<DNA_SIMPLE>;
    template class FastxParser<DNAN_SIMPLE>;
    template class SequenceReader<parsing::FastxParser<DNA_SIMPLE>>;
    template class SequenceReader<parsing::FastxParser<DNAN_SIMPLE>>;
    template class SplitPairedReader<parsing::FastxParser<DNA_SIMPLE>>;
    template class SplitPairedReader<parsing::FastxParser<DNAN_SIMPLE>>;

}
