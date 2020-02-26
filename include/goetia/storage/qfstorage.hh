/**
 * (c) Camille Scott, 2019
 * File   : qfstorage.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */
/* qfstorage.hh -- goetia-modified oxli storage
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
#include <string>

#include "goetia/meta.hh"
#include "goetia/storage/storage.hh"

struct quotient_filter;
typedef quotient_filter QF;

namespace goetia {
namespace storage {

/*
 * \class QFStorage
 *
 * \brief A Quotient Filter storage
 */
class QFStorage : public Storage<uint64_t>,
                  public Tagged<QFStorage> {
protected:
    std::shared_ptr<QF> cf;
    int _size;

public:
  
  using Storage<uint64_t>::value_type;

  QFStorage(int size);

  ~QFStorage();

  static std::shared_ptr<QFStorage> build(int size) {
      return std::make_shared<QFStorage>(size);
  }
  std::shared_ptr<QFStorage> clone() const;

  const bool insert(value_type khash);
  const count_t insert_and_query(value_type khash);

  // get the count for the given k-mer hash.
  const count_t query(value_type khash) const;

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

  double estimated_fp() {
      double fp = n_occupied() / get_tablesizes()[0];
      fp = pow(fp, n_tables());
      return fp;
  }

  static std::shared_ptr<QFStorage> deserialize(std::ifstream& in) {
    return {};
  }

  void serialize(std::ofstream& out) {}


};


template<>
struct is_probabilistic<QFStorage> { 
      static const bool value = true;
};



template<>
struct is_counting<QFStorage> {
    static const bool value = true;
};

}
}

#endif
