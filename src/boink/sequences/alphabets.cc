/**
 * (c) Camille Scott, 2019
 * File   : alphabets.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 07.01.2020
 */

#include "boink/sequences/alphabets.hh"

template <> const std::string boink::Alphabet<boink::DNA_SIMPLE>::SYMBOLS = "ACGT";
template <> const std::string boink::Alphabet<boink::DNA_SIMPLE>::COMPLEMENTS = "TGCA";

template <> const std::string boink::Alphabet<boink::DNAN_SIMPLE>::SYMBOLS = "ACGTN";
template <> const std::string boink::Alphabet<boink::DNAN_SIMPLE>::COMPLEMENTS = "TGCAN";

template <> const std::string boink::Alphabet<boink::IUPAC_NUCL>::SYMBOLS = "ATUGCYRSWKMBDHVN";
template <> const std::string boink::Alphabet<boink::IUPAC_NUCL>::COMPLEMENTS = "TAACGRYSWMKVHDBN";
