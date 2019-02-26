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
#include <utility>

#include "boink/boink.hh"
#include "boink/minimizers.hh"
#include "boink/hashing/alphabets.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/hashshifter.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/kmeriterator.hh"

#include "rollinghash/cyclichash.h"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#include "bbhash/BooPHF.h"
#pragma GCC diagnostic pop

namespace boink {
namespace hashing {


struct Unikmer {
    hash_t   hash;
    uint64_t partition;

    Unikmer(hash_t hash, uint64_t partition)
        : hash(hash), partition(partition)
    {}

    Unikmer(hash_t hash)
        : hash(hash), partition(ULLONG_MAX)
    {}

    Unikmer()
        : hash(ULLONG_MAX), partition(ULLONG_MAX)
    {}

    const bool is_valid() const {
        return partition != ULLONG_MAX;
    }

    friend bool operator>(const Unikmer& lhs, const Unikmer& rhs) {
        return lhs.hash > rhs.hash;
    }

    friend bool operator==(const Unikmer& lhs, const Unikmer& rhs) {
        return lhs.hash == rhs.hash;
    }

    friend std::ostream& operator<<(std::ostream& o, const Unikmer& un);
};


inline std::ostream& operator<<(std::ostream& o, const Unikmer& un) {
    o << "<Unikmer hash=" << un.hash;
    if (un.is_valid()) {
        o << " p=" << un.partition;
    } else {
        o << " p=None";
    }
    o << ">";
    return o;
}


class UKHS : public KmerClient {

protected:

    typedef boomphf::SingleHashFunctor<uint64_t>        boophf_hasher_t;
    typedef boomphf::mphf<uint64_t, boophf_hasher_t>    boophf_t;

    std::vector<uint64_t>          ukhs_hashes;
    std::vector<uint64_t>          ukhs_revmap;
    RollingHashShifter<DNA_SIMPLE> ukhs_hasher;
    bool                           ukhs_initialized;
    std::unique_ptr<boophf_t>      bphf;

public:

    explicit UKHS(uint16_t K,
                  std::vector<std::string>& ukhs)
        : KmerClient  (K),
          ukhs_revmap (ukhs.size()),
          ukhs_hasher (K)
    {
        if (ukhs.front().size() != K) {
            throw BoinkException("K does not match k-mer size from provided UKHS");
        }

        std::cerr << "Building MPHF with "
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

    bool query(Unikmer& unikmer) {
        unikmer.partition = ULLONG_MAX;
        unikmer.partition = bphf->lookup(unikmer.hash);
        if (!unikmer.is_valid()) {
            return false;
        }
        if (ukhs_revmap[unikmer.partition] == unikmer.hash) {
            return true;
        }
        return false;
    }

    std::vector<uint64_t> get_hashes() const {
        return ukhs_hashes;
    }

    const size_t n_hashes() const {
        return ukhs_hashes.size();
    }
};


template <const std::string& Alphabet = DNA_SIMPLE>
class UKHShifter : public HashShifter<UKHShifter<Alphabet>,
                                      Alphabet> {

protected:

    typedef HashShifter<UKHShifter<Alphabet>, Alphabet> BaseShifter;

    using BaseShifter::_K;

    CyclicHash<hash_t>           window_hasher;
    RollingHashShifter<Alphabet> ukhs_hasher;

    uint16_t                     window_K;
    uint16_t                     seed_K;

    InteriorMinimizer<Unikmer>   minimizer;
    std::shared_ptr<UKHS>        ukhs;

    bool                         ukhs_hasher_on_left;

    void update_unikmer() {
        Unikmer unikmer(ukhs_hasher.get());
        if (ukhs->query(unikmer)) {
            minimizer.update(unikmer);
        } else {
            minimizer.update(Unikmer());
        }
    }

public:

    typedef PartitionedHash hash_type;

    explicit UKHShifter(uint16_t K,
                        uint16_t seed_K,
                        std::shared_ptr<UKHS> ukhs)
        : BaseShifter   (K),
          window_hasher (K),
          ukhs_hasher   (seed_K),
          window_K      (K),
          seed_K        (seed_K),
          minimizer     (K - seed_K + 1),
          ukhs          (ukhs)
    {

    }

    explicit UKHShifter(UKHShifter& other)
        : UKHShifter(other.K(), other.ukhs_K(), other.get_ukhs())
    {
    }

    void init() {
        if (this->initialized) {
            return;
        }
        for (auto c : this->symbol_deque) {
            this->_validate(c);
            window_hasher.eat(c);
        }
        this->initialized = true;
    }

    const uint16_t ukhs_K() const {
        return seed_K;
    }

    hash_t get() {
        return window_hasher.hashvalue;
    }

    hash_t _hash(const std::string& sequence) const {
        return hash_cyclic(sequence, window_K);
    }

    hash_t _hash(const char * sequence) const {
        CyclicHash<hash_t> tmp_hasher(window_K);
        for (uint16_t i = 0; i < window_K; ++i) {
            tmp_hasher.eat(sequence[i]);
        }
        return tmp_hasher.hashvalue;
    }

    std::shared_ptr<UKHS> get_ukhs() {
        return ukhs;
    }

    std::vector<uint64_t> get_ukhs_hashes() const {
        return ukhs->get_hashes();
    }

    const size_t n_ukhs_hashes() const {
        return ukhs->n_hashes();
    }

    auto get_unikmers() const {
        return minimizer.get_minimizers();
    }

    auto get_front_unikmer() const {
        return minimizer.get_front_value();
    }

    auto get_front_partition() const {
        return minimizer.get_front_value().partition;
    }

    auto get_back_partition() const {
        return minimizer.get_back_value().partition;
    }
    
    bool query_unikmer(Unikmer& unikmer) {
        return ukhs->query(unikmer);
    }

    void reset_unikmers() {
        minimizer.reset();

        ukhs_hasher.set_cursor(this->symbol_deque.begin(),
                               this->symbol_deque.begin() + seed_K);
        
        for (auto cit = this->symbol_deque.begin() + seed_K;
                  cit != this->symbol_deque.end(); ++cit) {

            update_unikmer();
            ukhs_hasher.shift_right(*cit);
        }
        update_unikmer();

        if (minimizer.size() != 1) {
            throw BoinkException("Window should contain a k-mer from the UKHS!");
        }

        ukhs_hasher_on_left = false;

        init();
    }

    hash_t update_left(const char c) {
        // update main window left
        this->window_hasher.reverse_update(c, this->symbol_deque.back());

        // update ukhs hasher
        if (!ukhs_hasher_on_left) {
            // if we're not on the left of the window, reset the cursor there
            ukhs_hasher.set_cursor(this->symbol_deque.begin(), this->symbol_deque.begin() + seed_K);
            ukhs_hasher_on_left = true;
        } else {
            // othewise just shift the new symbol on
            ukhs_hasher.shift_left(c);
        }
        update_unikmer();
        return this->get();
    }

    hash_t update_right(const char c) {
        // update main window right
        this->window_hasher.update(this->symbol_deque.front(), c);

        // if the ukhs hasher is on the left, reset the cursor
        if (ukhs_hasher_on_left) {
            ukhs_hasher.set_cursor(this->symbol_deque.end() - seed_K, this->symbol_deque.end());
            ukhs_hasher_on_left = false;
        } else {
            ukhs_hasher.shift_right(c);
        }
        update_unikmer();
        return this->get();
    }

    std::vector<shift_t> gather_left() {
        std::vector<shift_t> hashes;
        const char back = this->symbol_deque.back();
        for (auto symbol : Alphabet) {
            window_hasher.reverse_update(symbol, back);
            shift_t result(window_hasher.hashvalue, symbol);
            hashes.push_back(result);
            window_hasher.update(symbol, back);
        }

        return hashes;
    }

    std::vector<shift_t> gather_right() {
        std::vector<shift_t> hashes;
        const char front = this->symbol_deque.front();
        for (auto symbol : Alphabet) {
            window_hasher.update(front, symbol);
            hashes.push_back(shift_t(window_hasher.hashvalue, symbol));
            window_hasher.reverse_update(front, symbol);
        }
        return hashes;
    }
};

typedef UKHShifter<DNA_SIMPLE> DefaultUKHSShifter;


// specializations for KmerIterator

template <>
template <>
typename hash_return<PartitionedHash>::type
KmerIterator<DefaultUKHSShifter>::first<PartitionedHash>();


template <>
template <>
typename hash_return<PartitionedHash>::type
KmerIterator<DefaultUKHSShifter>::next<PartitionedHash>();

}
}

#endif
