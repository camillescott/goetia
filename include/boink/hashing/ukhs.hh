/* ukhs.hh -- k-mer hash functions
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_UKHS_HH
#define BOINK_UKHS_HH

#include <climits>
#include <iostream>

#include "boink/boink.hh"
#include "boink/hashing/alphabets.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/hashshifter.hh"
#include "boink/hashing/rollinghashshifter.hh"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wall"
#include "bbhash/BooPHF.h"
#pragma GCC diagnostic pop

namespace boink {
namespace hashing {

template <const std::string& Alphabet = DNA_SIMPLE>
class UKHS : public RollingHashShifter<Alphabet> {

protected:

    typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
    typedef boomphf::mphf<uint64_t, hasher_t>    boophf_t;
    typedef RollingHashShifter<Alphabet>         BaseShifter;

    using BaseShifter::_K;

    uint16_t                     seed_K;
    std::vector<uint64_t>        ukhs_hashes;
    std::vector<uint64_t>        ukhs_revmap;
    RollingHashShifter<Alphabet> ukhs_hasher;
    bool                         ukhs_initialized;
    std::unique_ptr<boophf_t>    bphf;

    uint64_t                     cur_unikmer;
    uint64_t                     cur_partition;

    bool                         ukhs_hasher_on_left;

    void update_unikmer(uint64_t unikmer_hash) {
        uint64_t partition;
        std::cout << "Check " << unikmer_hash << " against " << cur_unikmer << std::endl;
        if (query(unikmer_hash, partition)) {
            std::cout << "Got " << partition << std::endl;
            if (unikmer_hash < cur_unikmer) {
                cur_unikmer = unikmer_hash;
                cur_partition = partition;
            }
        }
    }

public:

    explicit UKHS(uint16_t K,
                  uint16_t seed_K,
                  std::vector<std::string>& ukhs)
        : RollingHashShifter<Alphabet>(K),
          seed_K                      (seed_K),
          ukhs_revmap                 (ukhs.size()),
          ukhs_hasher                 (seed_K),
          ukhs_initialized            (false),
          cur_unikmer                 (ULLONG_MAX),
          cur_partition               (ULLONG_MAX)
    {
        if (ukhs.front().size() != seed_K) {
            throw BoinkException("K does not match k-mer size from provided UKHS");
        }

        std::cerr << "Building MPHF for (W=" << K << ",K=" 
                  << seed_K << ") UKHS with "
                  << ukhs.size() << " k-mers."
                  << std::endl;

        for (auto unikmer : ukhs) {
            ukhs_hashes.push_back(ukhs_hasher.hash(unikmer));
        }
        bphf = std::make_unique<boophf_t>(ukhs_hashes.size(), ukhs_hashes, 1, 2.0);
        for (auto unikmer_hash : ukhs_hashes) {
            ukhs_revmap[bphf->lookup(unikmer_hash)] = unikmer_hash;
        }
        std::cerr << "Finished building MPHF." << std::endl;
    }

    std::vector<uint64_t> get_hashes() const {
        return ukhs_hashes;
    }

    bool query(uint64_t unikmer_hash, uint64_t& partition) {
        partition = ULLONG_MAX;
        partition = bphf->lookup(unikmer_hash);
        if (partition == ULLONG_MAX) {
            return false;
        }
        if (ukhs_revmap[partition] == unikmer_hash) {
            return true;
        }
        return false;
    }

    hash_t get_unikmer() const {
        return cur_unikmer;
    }

    hash_t get_partition() const {
        return cur_partition;
    }

    void reset_unikmer() {
        ukhs_hasher.set_cursor(this->symbol_deque.begin(),
                               this->symbol_deque.begin() + seed_K);
        std::cout << ukhs_hasher.get_cursor() << std::endl;
        for (auto cit = this->symbol_deque.begin() + seed_K;
                  cit != this->symbol_deque.end(); ++cit) {

            update_unikmer(ukhs_hasher.get());
            ukhs_hasher.shift_right(*cit);
        }
        update_unikmer(ukhs_hasher.get());

        if (cur_unikmer == ULLONG_MAX) {
            throw BoinkException("Window should contain a k-mer from the UKHS!");
        }

        ukhs_hasher_on_left = false;

        RollingHashShifter<Alphabet>::init();
    }

    hash_t update_left(const char c) {
        // update main window left
        this->hasher.reverse_update(c, this->symbol_deque.back());

        // update ukhs hasher
        if (!ukhs_hasher_on_left) {
            // if we're not on the left of the window, reset the cursor there
            ukhs_hasher.set_cursor(this->symbol_deque.begin(), this->symbol_deque.begin() + seed_K);
            ukhs_hasher_on_left = true;
        } else {
            // othewise just shift the new symbol on
            ukhs_hasher.shift_left(c);
        }
        update_unikmer(ukhs_hasher.hashvalue);
        return this->get();
    }

    hash_t update_right(const char c) {
        // update main window right
        this->hasher.update(this->symbol_deque.front(), c);

        // if the ukhs hasher is on the left, reset the cursor
        if (ukhs_hasher_on_left) {
            ukhs_hasher.set_cursor(this->symbol_deque.end() - seed_K, this->symbol_deque.end());
            ukhs_hasher_on_left = false;
        } else {
            ukhs_hasher.shift_right(c);
        }
        update_unikmer(ukhs_hasher.hashvalue);
        return this->get();
    }
};


typedef UKHS<DNA_SIMPLE> UKHSShifter;


}
}

#endif
