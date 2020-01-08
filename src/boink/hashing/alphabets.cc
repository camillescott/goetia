/**
 * (c) Camille Scott, 2019
 * File   : alphabets.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 06.08.2019
 */
/* alphabets.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 * 
 * This file is part of khmer, https://github.com/dib-lab/khmer/, and is
 * Copyright (C) 2015-2016, The Regents of the University of California.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of the Michigan State University nor the names
 *       of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * LICENSE (END)
 *
 * Contact: khmer-project@idyll.org
 */

#include "boink/hashing/alphabets.hh"
#include "boink/hashing/rollinghash/cyclichash.h"

#define tbl \
  "                                                                "\
  /*ABCDEFGHIJKLMNOPQRSTUVWXYZ      abcdefghijklmnopqrstuvwxyz    */\
  " TVGH FCD  M KN   YSAABW R       TVGH FCD  M KN   YSAABW R"
  //" TVGH FCD  M KA   YSAABWARA      TVGH FCD  M KA   YSAABWARA"

namespace boink {
namespace hashing {

std::string DNA_SIMPLE = "ACGT";
std::string DNA_SIMPLE_CMP = "TGCA";

std::string DNAN_SIMPLE = "ACGTN";
std::string DNAN_SIMPLE_CMP = "TGCAN";

std::string RNA_SIMPLE = "ACGUT";
std::string RNAN_SIMPLE = "ACGUTN";

std::string IUPAC_NUCL = "ACGTURYSWKMBDHVN.-";
std::string IUPAC_AA = "ACDEFGHIKLMNPQRSTVWY";


std::string revcomp(const std::string& kmer) {
    std::string out = kmer;

    auto from = out.begin();
    auto to = out.end();

    char c;
    for (to--; from <= to; from++, to--) {
        c = tbl[(int)*from];
        *from = tbl[(int)*to];
        *to = c;
    }
    return out;
}


uint64_t hash_cyclic(const std::string& kmer, const uint16_t k) {
    CyclicHash<uint64_t> hasher(k);
    for (uint16_t i = 0; i < k; ++i) {
        hasher.eat(kmer[i]);
    }
    return hasher.hashvalue;
}


}
}
