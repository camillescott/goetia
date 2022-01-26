/**
(c) Camille Scott, 2021
File   : hllcounter.hh
License: MIT
Author : Camille Scott <camille.scott.w@gmail.com>
Date   : 15.03.2021

This file modified on the original found in khmer.

---

This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2014-2015, Michigan State University.
Copyright (C) 2015-2016, The Regents of the University of California.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.
    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)
Contact: khmer-project@idyll.org
*/

#ifndef HLLCOUNTER_HH
#define HLLCOUNTER_HH

#include <cstdint>
#include <string>
#include <vector>

#include "goetia/goetia.hh"
#include "goetia/meta.hh"
#include "goetia/storage/storage.hh"
#include "goetia/sketches/sketch/include/sketch/hll.h"

namespace goetia {

class HLLStorage;

template<>
struct StorageTraits<HLLStorage> {
    static constexpr bool is_probabilistic = true;
    static constexpr bool is_counting      = true;

    typedef std::tuple<double> params_type;
    static constexpr params_type default_params = std::make_tuple(0.05);
};


class HLLStorage : public Storage<uint64_t>,
                   public Tagged<HLLStorage>
{

/*
 * This is just a thin wrapper around sketch::hll_t to make it 
 * compatible with partitioned_storage for use as a backend
 * in UnikmerSketch. It is *not* meant to be used as an
 * actual dBG storage.
 */
public:

    typedef sketch::hll::shll_t sketch_t;

protected:

    std::shared_ptr<sketch_t> counter;

public:

    using Storage<uint64_t>::value_type;
    using Traits = StorageTraits<HLLStorage>;
    const double error_rate;

    HLLStorage(double error_rate)
        : error_rate(error_rate) {
        int p = ceil(log2(pow(1.04 / error_rate, 2)));
        counter = std::make_shared<sketch_t>(p);
    }

    //HLLStorage(int nc) {
    //    counter = std::make_shared<sketches::HLLCounter>(nc);
    //}

    std::shared_ptr<HLLStorage> clone() const {
        return std::make_shared<HLLStorage>(error_rate);
    }

    static std::shared_ptr<HLLStorage> build(double error_rate) {
        return std::make_shared<HLLStorage>(error_rate);
    }

    static std::shared_ptr<HLLStorage> build(const typename Traits::params_type& params) {
        return make_shared_from_tuple<HLLStorage>(params);
    }

    void reset() {
        counter = std::make_shared<sketch_t>(error_rate);
    }

    const uint64_t n_occupied() const {
        return 0;
    }

    const uint64_t n_unique_kmers() const {
        return counter->report();
    }

    const inline bool insert(value_type h) {
        counter->addh(h);
        return true;
    }

    const count_t insert_and_query(value_type h) {
        counter->addh(h);
        return 1;
    }

    const count_t query(value_type h) const {
        return 1;
    }

    void save(std::string, uint16_t ) {
    }

    void load(std::string, uint16_t &) {
    }

    static std::shared_ptr<HLLStorage> deserialize(std::ifstream& in) {
    }

    void serialize(std::ofstream& out) {
    }

    byte_t ** get_raw_tables() {
        return nullptr;
    }
};


}

#endif // HLLCOUNTER_HH

