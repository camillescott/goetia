/**
 * (c) Camille Scott, 2019
 * File   : utagger.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 05.08.2019
 */

#ifndef BOINK_UTAGGER_HH
#define BOINK_UTAGGER_HH


#include "boink/traversal.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/dbg.hh"
#include "boink/cdbg/udbg.hh"

namespace boink {
namespace cdbg {

template <class GraphType>
struct UTagger {
    using ShifterType     = typename GraphType::shifter_type;
    using TraversalType   = Traverse<GraphType>;
    using TraverserType   = typename TraversalType::dBG;
    using MinimizerType   = WKMinimizer<ShifterType>;
    using UnikmerIterType = typename hashing::KmerIterator<hashing::UKHShifter>;

    using cDBGType        = uDBG<GraphType>;
    using TagType         = typename cDBGType::Tag;


};
}
}

#endif
