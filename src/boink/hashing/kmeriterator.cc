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
#include "boink/hashing/ukhs.hh"
#include "boink/storage/storage_types.hh"
#include "boink/dbg.hh"
#include "boink/traversal.hh"

namespace boink {
namespace hashing {

template<>
template<>
KmerIterator<UKHS::LazyShifter>::KmerIterator(const std::string& seq, uint16_t K)
    : KmerClient(K) {
    // this is dumb and should be compile time but i don't feel like figuring it out
    throw BoinkException("Invalid constructor for LazyShifter (need W and K");
}


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

    auto ret = shifter->shift_right(_seq[index + _K - 1]);
    index += 1;

    return ret;
}


template<class ShifterType>
typename ShifterType::hash_type
KmerIterator<ShifterType>::first() {
    _initialized = true;

    index += 1;
    return shifter->set_cursor(_seq);
}


template<class ShifterType>
bool
KmerIterator<ShifterType>::done() const {
    return (index + _K > _seq.length());
}



}
}

template class boink::hashing::KmerIterator<boink::hashing::RollingHashShifter>;
template class boink::hashing::KmerIterator<boink::hashing::UKHS::LazyShifter>;

template class boink::hashing::KmerIterator<boink::Traverse<boink::dBG<boink::storage::SparseppSetStorage,
                                                                       boink::hashing::RollingHashShifter>>::dBG>;
template class boink::hashing::KmerIterator<boink::Traverse<boink::dBG<boink::storage::BitStorage,
                                                                       boink::hashing::RollingHashShifter>>::dBG>;
template class boink::hashing::KmerIterator<boink::Traverse<boink::dBG<boink::storage::ByteStorage,
                                                                       boink::hashing::RollingHashShifter>>::dBG>;
template class boink::hashing::KmerIterator<boink::Traverse<boink::dBG<boink::storage::NibbleStorage,
                                                                       boink::hashing::RollingHashShifter>>::dBG>;
template class boink::hashing::KmerIterator<boink::Traverse<boink::dBG<boink::storage::QFStorage,
                                                                       boink::hashing::RollingHashShifter>>::dBG>;


template class boink::hashing::KmerIterator<boink::Traverse<boink::dBG<boink::storage::SparseppSetStorage,
                                                                       boink::hashing::UKHS::LazyShifter>>::dBG>;
template class boink::hashing::KmerIterator<boink::Traverse<boink::dBG<boink::storage::BitStorage,
                                                                       boink::hashing::UKHS::LazyShifter>>::dBG>;
template class boink::hashing::KmerIterator<boink::Traverse<boink::dBG<boink::storage::ByteStorage,
                                                                       boink::hashing::UKHS::LazyShifter>>::dBG>;
template class boink::hashing::KmerIterator<boink::Traverse<boink::dBG<boink::storage::NibbleStorage,
                                                                       boink::hashing::UKHS::LazyShifter>>::dBG>;
template class boink::hashing::KmerIterator<boink::Traverse<boink::dBG<boink::storage::QFStorage,
                                                                       boink::hashing::UKHS::LazyShifter>>::dBG>;
