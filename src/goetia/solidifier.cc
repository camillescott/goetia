/**
 * (c) Camille Scott, 2019
 * File   : solidifier.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 12.03.2020
 */

#include "goetia/solidifier.hh"

namespace goetia {

template class StreamingSolidFilter<dBG<storage::BitStorage, hashing::FwdLemireShifter>>;
template class StreamingSolidFilter<dBG<storage::QFStorage, hashing::FwdLemireShifter>>;

}
