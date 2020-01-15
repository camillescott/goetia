/**
 * (c) Camille Scott, 2019
 * File   : ukhs.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 06.08.2019
 */


#ifndef BOINK_UKHS_HH
#define BOINK_UKHS_HH

#include <algorithm>
#include <climits>
#include <iostream>
#include <memory>
#include <optional>
#include <utility>

#include "boink/boink.hh"
#include "boink/meta.hh"
#include "boink/sequences/alphabets.hh"
#include "boink/hashing/canonical.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/kmers/kmerclient.hh"
#include "boink/storage/sparsepp/spp.h"

#include "rollinghash/cyclichash.h"


namespace boink::hashing { 

template <class ShifterType = RollingHashShifter<>>
struct UKHS : public kmers::KmerClient {

    typedef typename ShifterType::hash_type hash_type;
    typedef typename hash_type::value_type  value_type;

    typedef Partitioned<hash_type>          Unikmer;


protected:

    // map of unikmer hashes to partitions
    typedef spp::sparse_hash_map<value_type, uint64_t> pmap_t;
    pmap_t                                             pmap;
    std::vector<hash_type>                             hashes;


public:

    // Window size W. K size is stored in KmerClient.
    const uint16_t _W;

    explicit UKHS(uint16_t W,
                 uint16_t K,
                 std::vector<std::string>& ukhs);

    const uint16_t W() const {
        return _W;
    }

    static std::shared_ptr<UKHS> build(uint16_t W,
                                      uint16_t K,
                                      std::vector<std::string>& ukhs) {
        return std::make_shared<UKHS>(W, K, ukhs);
    }

    std::optional<Unikmer> query(hash_type unikmer_hash);

    std::vector<hash_type> get_hashes() const {
        return hashes;
    }

    const size_t n_hashes() const {
        return pmap.size();
    }
};

}

#endif
