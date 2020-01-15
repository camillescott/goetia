/**
 * (c) Camille Scott, 2019
 * File   : kmeriterator.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 29.08.2019
 */


#include <string>

#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/ukhshashshifter.hh"
#include "boink/hashing/hashextender.hh"
#include "boink/storage/storage_types.hh"
#include "boink/traversal.hh"
#include "boink/dbg.hh"

namespace boink::hashing {

    /*
template<>
template<>
KmerIterator<UKHS::LazyShifter>::KmerIterator(const std::string& seq, uint16_t K)
    : KmerClient(K) {
    // this is dumb and should be compile time but i don't feel like figuring it out
    throw BoinkException("Invalid constructor for LazyShifter (need W and K");
}
*/

template<class ShifterType>
KmerIterator<ShifterType>::KmerIterator(const std::string& seq, ShifterType * shifter)
    : KmerClient(shifter->K()), 
      _seq(seq),
      index(0), 
      _initialized(false),
      _shifter_owner(false), 
      shifter(shifter) 
{
    if (_seq.length() < _K) {
        throw SequenceLengthException("Sequence must have length >= K");
    }
}


template<class ShifterType>
KmerIterator<ShifterType>::KmerIterator(const std::string& seq, ShifterType& shifter_proto)
    : KmerClient(shifter_proto.K()), 
      _seq(seq), 
      index(0), 
      _initialized(false), 
      _shifter_owner(true)
{

    if (_seq.length() < _K) {
        throw SequenceLengthException("Sequence must have length >= K");
    }
    shifter = new ShifterType(shifter_proto);
}



template<class ShifterType>
typename ShifterType::hash_type
KmerIterator<ShifterType>::next() {
    if (!_initialized) {
        return first();
    }

    if (done()) {
        throw InvalidCharacterException("past end of iterator");
    }

    auto ret = shifter->shift_right(_seq[index - 1], _seq[index + _K - 1]);
    index += 1;

    return ret;
}


template<class ShifterType>
typename ShifterType::hash_type
KmerIterator<ShifterType>::first() {
    _initialized = true;

    index += 1;
    return shifter->hash_base(_seq);
}


template<class ShifterType>
bool
KmerIterator<ShifterType>::done() const {
    return (index + _K > _seq.length());
}

template class KmerIterator<FwdRollingShifter>;
template class KmerIterator<CanRollingShifter>;

template class KmerIterator<FwdUnikmerShifter>;
template class KmerIterator<CanUnikmerShifter>;

template class KmerIterator<HashExtender<FwdRollingShifter>>;
template class KmerIterator<HashExtender<CanRollingShifter>>;

template class KmerIterator<HashExtender<FwdUnikmerShifter>>;
template class KmerIterator<HashExtender<CanUnikmerShifter>>;

}

namespace boink {

template class hashing::KmerIterator<dBGWalker<dBG<storage::BitStorage, hashing::FwdRollingShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::BitStorage, hashing::CanRollingShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::BitStorage, hashing::FwdUnikmerShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::BitStorage, hashing::CanUnikmerShifter>>>;

template class hashing::KmerIterator<dBGWalker<dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::SparseppSetStorage, hashing::CanRollingShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>>>;

template class hashing::KmerIterator<dBGWalker<dBG<storage::ByteStorage, hashing::FwdRollingShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::ByteStorage, hashing::CanRollingShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::ByteStorage, hashing::FwdUnikmerShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::ByteStorage, hashing::CanUnikmerShifter>>>;

template class hashing::KmerIterator<dBGWalker<dBG<storage::NibbleStorage, hashing::FwdRollingShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::NibbleStorage, hashing::CanRollingShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::NibbleStorage, hashing::CanUnikmerShifter>>>;

template class hashing::KmerIterator<dBGWalker<dBG<storage::QFStorage, hashing::FwdRollingShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::QFStorage, hashing::CanRollingShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::QFStorage, hashing::FwdUnikmerShifter>>>;
template class hashing::KmerIterator<dBGWalker<dBG<storage::QFStorage, hashing::CanUnikmerShifter>>>;

}

