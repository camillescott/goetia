/* qfstorage.hh -- boink-modified oxli storage
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 *
 *** END BOINK LICENSE BLOCK
 *
 * This file is part of khmer, https://github.com/dib-lab/khmer/, and is
 * Copyright (C) 2016, The Regents of the University of California.
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
 *     * Neither the name of the University of California nor the names
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

#ifndef BOINK_QFSTORAGE_HH
#define BOINK_QFSTORAGE_HH

#include <cassert>
#include <array>
#include <memory>
#include <mutex>
#include <unordered_map>

#include "boink/hashing/hashing_types.hh"
#include "boink/storage/storage.hh"


namespace boink {
namespace storage {

typedef struct quotient_filter;
typedef quotient_filter QF;

/*
 * \class QFStorage
 *
 * \brief A Quotient Filter storage
 */
 class QFStorage : public Storage {
protected:
    std::shared_ptr<QF> cf;

public:
  QFStorage(int size);

  ~QFStorage();

  count_t test_and_set_bits(hashing::hash_t khash) {
    count_t x = get_count(khash);
    add(khash);
    return !x;
  }

  //
  bool add(hashing::hash_t khash);

  // get the count for the given k-mer hash.
  const count_t get_count(hashing::hash_t khash) const;

  // Accessors for protected/private table info members
  // xnslots is larger than nslots. It includes some extra slots to deal
  // with some details of how the counting is implemented
  std::vector<uint64_t> get_tablesizes() const;
  const size_t n_tables() const { return 1; }
  const uint64_t n_unique_kmers() const;
  const uint64_t n_occupied() const;
  void save(std::string outfilename, uint16_t ksize);
  void load(std::string infilename, uint16_t &ksize);

  byte_t **get_raw_tables() { return nullptr; }
  void reset() {}; //nop
};

}
}

#endif
