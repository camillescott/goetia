/**
 * (c) Camille Scott, 2019
 * File   : ukhs.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 06.08.2019
 */
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

#include <algorithm>
#include <climits>
#include <iostream>
#include <memory>
#include <optional>
#include <utility>

#include "boink/boink.hh"
#include "boink/meta.hh"
#include "boink/hashing/alphabets.hh"
#include "boink/hashing/hashshifter.hh"
#include "boink/hashing/rollinghash/cyclichash.h"
#include "boink/hashing/kmeriterator.hh"
#include "boink/kmers/kmerclient.hh"

#include "rollinghash/cyclichash.h"


namespace boomphf {

template <typename Item>
class SingleHashFunctor;

template <typename elem_t, typename Hasher_t>
class mphf;

}

namespace boink {
namespace hashing {

struct UKHS {

    using value_type = uint64_t;

    struct Unikmer {
        /* A hash and its associated partition. 
         *
         */
        value_type hash;
        uint64_t   partition;

        Unikmer(value_type hash, uint64_t partition)
            : hash(hash), partition(partition)
        {}

        Unikmer(value_type hash)
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

        friend bool operator<(const Unikmer& lhs, const Unikmer& rhs) {
            return lhs.hash < rhs.hash;
        }

        friend bool operator==(const Unikmer& lhs, const Unikmer& rhs) {
            return lhs.hash == rhs.hash;
        }

        friend inline std::ostream&
        operator<<(std::ostream& o, const Unikmer& un) {
            o << "<Unikmer hash=" << un.hash;
            if (un.is_valid()) {
                o << " p=" << un.partition;
            } else {
                o << " p=None";
            }
            o << ">";
            return o;
        }
    };

    struct BinnedKmer {
        value_type hash;
        Unikmer    unikmer;

        BinnedKmer(value_type hash, Unikmer unikmer)
            : hash(hash), unikmer(unikmer)
        {
        }

        BinnedKmer(value_type hash)
            : hash(hash)
        {
        }

        BinnedKmer()
            : hash(ULLONG_MAX)
        {
        }

        BinnedKmer(BinnedKmer& other)
            : hash(other.hash),
              unikmer(other.unikmer)
        {
        }

        BinnedKmer(const BinnedKmer& other)
            : hash(other.hash),
              unikmer(other.unikmer)
        {
        }

        operator value_type() const {
            return hash;
        }

        friend bool operator>(const BinnedKmer& lhs, const BinnedKmer& rhs) {
            return lhs.hash > rhs.hash;
        }

        friend bool operator<(const BinnedKmer& lhs, const BinnedKmer& rhs) {
            return lhs.hash < rhs.hash;
        }

        friend bool operator==(const BinnedKmer& lhs, const BinnedKmer& rhs) {
            return lhs.hash == rhs.hash;
        }

        friend inline std::ostream&
        operator<<(std::ostream& o, const BinnedKmer& kmer) {
            o << "<BinnedKmer hash=" << kmer.hash
              << " unikmer="         << kmer.unikmer
              << ">";
            return o;
        }
    };

    class Map : public kmers::KmerClient {

    protected:

        typedef boomphf::SingleHashFunctor<value_type>     boophf_hasher_t;
        typedef boomphf::mphf<value_type, boophf_hasher_t> boophf_t;

        std::vector<value_type>   ukhs_hashes;
        std::vector<value_type>   ukhs_revmap;
        bool                      ukhs_initialized;
        std::unique_ptr<boophf_t> bphf;

    public:

        // Window size W. K size is stored in KmerClient.
        const uint16_t _W;

        explicit Map(uint16_t W,
                     uint16_t K,
                     std::vector<std::string>& ukhs);

        ~Map();

        const uint16_t W() const {
            return _W;
        }

        static std::shared_ptr<Map> build(uint16_t W,
                                          uint16_t K,
                                          std::vector<std::string>& ukhs) {
            return std::make_shared<Map>(W, K, ukhs);
        }

        value_type query_revmap(uint64_t partition) {
            if (partition > ukhs_revmap.size()) {
                return ULLONG_MAX;
            }
            return ukhs_revmap[partition];
        }

        bool query(Unikmer& unikmer);

        std::vector<value_type> get_hashes() const {
            return ukhs_hashes;
        }

        const size_t n_hashes() const {
            return ukhs_hashes.size();
        }
    };

    class LazyShifter : public HashShifter<LazyShifter, BinnedKmer>,
                        public Tagged<LazyShifter> {

        /* Shifter that keeps track of each k-mer's associated
         * Unikmer. The hash_type will end being composed as:
         *
         *     hash_type:
         *         uint64_t hash
         *         Unikmer  unikmer:
         *             uint64_t hash
         *             size_t   partition
         *
         * There will usually be only one associated unikmer
         * for a single k-mer. Sometimes, however, a few unikmers
         * will exist within the window; in this case, we choose
         * the minimum unikmer as the representative. Rather than
         * use the minizer class for this, we'll just keep two deques,
         * one with the indices and one with the unikmers, to track
         * the unikmers within the window; the logic is easier 
         * this way concerning directionality, and there will usually
         * only be one unikmer in the window anyway, so it will have
         * less overhead.
         */

        typedef HashShifter<LazyShifter, BinnedKmer> BaseShifter;

    public:

        typedef BaseShifter baseshifter_type;

        // These types will now be composited with Unikmer
        using BaseShifter::hash_type;
        using BaseShifter::shift_type;
        using BaseShifter::kmer_type;

    protected:

        using BaseShifter::_K;
        using BaseShifter::kmer_window;

        // K size of the window, "W" of the W,K UKHS
        // Same as HashShifter->_K
        uint16_t _window_K;
        // K size of the unikmer
        uint16_t _unikmer_K;

        // Unikmer uses value_type from rolling hashshifter
        CyclicHash<value_type>       window_hasher;
        CyclicHash<value_type>       unikmer_hasher;

        // Current position of the unikmer hasher within
        // the window: it will need to be moved across the window
        // if we change directions
        bool unikmer_hasher_on_left;

        // We'll store the unikmers for the current window and their
        // indices, only finding the minimum when queried.
        std::deque<Unikmer>    window_unikmers;
        std::deque<size_t>     unikmer_indices;

        // Reset the unikmer hasher and eat the length K-1 window prefix
        void eat_unikmer_hasher_left_prefix() {
            for (uint16_t i = 0; i < _unikmer_K - 1; ++i) {
                unikmer_hasher.eat(*(kmer_window.begin() + i));
            }
        }

        // Reset the unikmer hasher and eat the length K-1 window suffix
        void eat_unikmer_hasher_right_suffix() {
            for (uint16_t i = this->_K - _unikmer_K + 1; i < this->_K; ++i) {
                unikmer_hasher.eat(*(kmer_window.begin() + i));
            }
        }

        void shift_unikmers_left() {
            for (auto& index : unikmer_indices) {
                ++index;
            }
            if (unikmer_indices.back() > this->_K - _unikmer_K) {
                unikmer_indices.pop_back();
                window_unikmers.pop_back();
            }
        }

        void update_unikmer_left(const char c) {
            if (!unikmer_hasher_on_left) {
                // if we're not on the left of the window, reset the cursor there
                unikmer_hasher.reset();
                unikmer_hasher.eat(c);
                eat_unikmer_hasher_left_prefix();
                unikmer_hasher_on_left = true;
            } else {
                // othewise just shift the new symbol on
                   //std::cout << "Reverse update, "
                //          << *(kmer_window.begin() + _unikmer_K - 1)
                //          << " => "
                //          << c << std::endl;
                unikmer_hasher.reverse_update(c, *(kmer_window.begin() + _unikmer_K - 1));
            }

            shift_unikmers_left();
            Unikmer cur_unikmer(unikmer_hasher.hashvalue);
            if (ukhs_map->query(cur_unikmer)) {
                unikmer_indices.push_front(0);
                window_unikmers.push_front(cur_unikmer);
            }
        }

        void shift_unikmers_right() {
            //std::cout << "indices: " << repr(unikmer_indices) << std::endl;
            //std::cout << "unikmers: " << repr(window_unikmers) << std::endl;
            if (unikmer_indices.front() == 0) {
                unikmer_indices.pop_front();
                window_unikmers.pop_front();
            }

            for (auto& index : unikmer_indices) {
                --index;
            }

            //std::cout << "new indices: " << repr(unikmer_indices) << std::endl;
            //std::cout << "new unikmers: " << repr(window_unikmers) << std::endl;
        }

        void update_unikmer_right(const char c) {
            //std::cout << "LazyShifter update_unikmer_right(" 
            //          << c << "): " << this->get_cursor() << std::endl;
            // if the ukhs hasher is on the left, reset the cursor
            if (unikmer_hasher_on_left) {
                //std::cout << "LazyShifter set_unikmer_hasher_right_suffix" << std::endl;
                unikmer_hasher.reset();
                eat_unikmer_hasher_right_suffix();
                unikmer_hasher.eat(c);
                unikmer_hasher_on_left = false;
            } else {
                //std::cout << "LazyShifter update: "
                //          << *(kmer_window.begin() + this->_K - _unikmer_K)
                //          << " => " << c << std::endl;
                unikmer_hasher.update(*(kmer_window.begin() + this->_K - _unikmer_K), c);
            }
            
            shift_unikmers_right();
            Unikmer cur_unikmer(unikmer_hasher.hashvalue);
            if (ukhs_map->query(cur_unikmer)) {
                //std::cout << "LazyShifter: found unikmer" << std::endl;
                unikmer_indices.push_back(this->_K - _unikmer_K);
                window_unikmers.push_back(cur_unikmer);
            }
            //std::cout << "LazyShifter: leave update_unikmer, uhash: " << cur_unikmer << std::endl;
        }

        Unikmer get_min_unikmer() {
            //std::cout << "LazyShifter: " << repr(window_unikmers) << std::endl;
            //std::cout << "LazyShifter: " << repr(unikmer_indices) << std::endl;

            if (window_unikmers.size() == 0) {
                throw BoinkException("Window should contain unikmer.");
            }
            return *std::min_element(std::begin(window_unikmers),
                                     std::end(window_unikmers));
        }

        void clear_unikmers() {
            unikmer_indices.clear();
            window_unikmers.clear();
        }

    public:

        using BaseShifter::symbols;

        // W,K unikmer map
        std::shared_ptr<Map>         ukhs_map;

        explicit LazyShifter(uint16_t K,
                             uint16_t _unikmer_K,
                             std::shared_ptr<Map> ukhs)
            : BaseShifter    (K),
              window_hasher  (K),
              unikmer_hasher (_unikmer_K),
              _window_K       (K),
              _unikmer_K      (_unikmer_K),
              ukhs_map       (ukhs)
        {
            if (ukhs_map->W() != K) {
                throw BoinkException("Shifter K does not match UKHS::Map W.");
            }
            if (ukhs_map->K() != _unikmer_K) {
                throw BoinkException("Shifter _unikmer_K does not match UKHS::Map K.");
            }
        }

        explicit LazyShifter(LazyShifter& other)
            : LazyShifter(other.window_K(),
                          other.unikmer_K(),
                          other.ukhs_map)
        {
        }

        explicit LazyShifter(const LazyShifter& other)
            : LazyShifter(other.window_K(),
                          other.unikmer_K(),
                          other.ukhs_map)
        {
        }

        uint16_t window_K() const {
            return _window_K;
        }

        uint16_t unikmer_K() const {
            return _unikmer_K;
        }


        void init() {
            if (this->initialized) {
                return;
            }
            //std::cout << "LazyShifter: init()" << std::endl;
            for (uint16_t i = 0; i < _unikmer_K; ++i) {
                this->_validate(*(kmer_window.begin() + i));
                window_hasher.eat(*(kmer_window.begin() + i));
                unikmer_hasher.eat(*(kmer_window.begin() + i));
                //std::cout << "LazyShifter: eat " << *(kmer_window.begin() +i) << std::endl;
            }

            Unikmer unikmer(unikmer_hasher.hashvalue);
            if (ukhs_map->query(unikmer)) {
                unikmer_indices.push_back(0);
                window_unikmers.push_back(unikmer);
            }

            for (uint16_t i = _unikmer_K; i < _window_K; ++i) {
                // once we've eaten the first K bases, start keeping
                // track of unikmers

                window_hasher.eat(*(kmer_window.begin() + i));
                unikmer_hasher.update(*(kmer_window.begin() + i - _unikmer_K),
                                      *(kmer_window.begin() + i));

                Unikmer unikmer(unikmer_hasher.hashvalue);
                if (ukhs_map->query(unikmer)) {
                    unikmer_indices.push_back(i - _unikmer_K + 1);
                    window_unikmers.push_back(unikmer);
                }
                //std::cout << "LazyShifter: eat/update "
                //          << *(kmer_window.begin() +i)
                //          << ", shiftoff " <<  *(kmer_window.begin() + i - _unikmer_K)
                //          << std::endl;
            }

            this->initialized = true;
            unikmer_hasher_on_left = false;
            //std::cout << "LazyShifter: finish init(): " << this->get() << std::endl;
        }

        std::vector<uint64_t> get_ukhs_hashes() const {
            return ukhs_map->get_hashes();
        }

        const size_t n_ukhs_hashes() const {
            return ukhs_map->n_hashes();
        }

        hash_type get() {
            return BinnedKmer(window_hasher.hashvalue, get_min_unikmer());
        }


        hash_type _hash(const std::string& sequence) const {
            LazyShifter shifter(this->_window_K,
                                this->_unikmer_K,
                                this->ukhs_map);
            shifter.set_cursor(sequence);
            return shifter.get();
        }

        hash_type _hash(const char * sequence) const {
            LazyShifter shifter(this->_window_K,
                                this->_unikmer_K,
                                this->ukhs_map);
            shifter.set_cursor(sequence);
            return shifter.get();
        }

        hash_type update_left(const char c) {
            // update main window left
            this->window_hasher.reverse_update(c, this->kmer_window.back());
            update_unikmer_left(c);
            return this->get();
        }

        hash_type update_right(const char c) {
            // update main window right
            this->window_hasher.update(this->kmer_window.front(), c);
            update_unikmer_right(c);
            return this->get();
        }

        std::vector<shift_type> gather_left() {
            std::vector<shift_type> hashes;

            // First get the min unikmer in the W-1 prefix, if there is one
            auto _last = unikmer_indices.back() > this->_K - _unikmer_K ? 
                         std::end(window_unikmers) - 1 :
                         std::end(window_unikmers);
            auto _min = std::min_element(std::begin(window_unikmers),
                                         _last);
            std::optional<Unikmer> current_min;
            if (_min != _last) {
                current_min = *_min;
            }
            
            if (!unikmer_hasher_on_left) {
                unikmer_hasher.reset();
                for (uint16_t i = 0; i < _unikmer_K; ++i) {
                    unikmer_hasher.eat(*(kmer_window.begin() + i));
                }
                unikmer_hasher_on_left = true;
            }

            // now compare each neighbor unikmer to the current min
            // to see if the neighbor w-mer has a different minimizer
            const char back = this->kmer_window.back();
            const char uback = *(this->kmer_window.begin() + _unikmer_K - 1);
            for (auto& symbol : symbols) {
                window_hasher.reverse_update(symbol, back);
                unikmer_hasher.reverse_update(symbol, uback);

                Unikmer cur_unikmer(unikmer_hasher.hashvalue);
                if (!current_min || (ukhs_map->query(cur_unikmer) && cur_unikmer < current_min.value())) {
                    hashes.push_back(shift_type(BinnedKmer(window_hasher.hashvalue, cur_unikmer),
                                     symbol));
                } else {
                    hashes.push_back(shift_type(BinnedKmer(window_hasher.hashvalue, current_min.value()),
                                     symbol));
                }

                // update hashers back to root 0
                window_hasher.update(symbol, back);
                unikmer_hasher.update(symbol, uback);

            }

            return hashes;
        }

        std::vector<shift_type> gather_right() {
            std::vector<shift_type> hashes;

            // First get the min unikmer in the W-1 prefix, if there is one
            auto _first = unikmer_indices.front() == 0 ? 
                          std::begin(window_unikmers) + 1 :
                          std::begin(window_unikmers);
            auto _min = std::min_element(_first, std::end(window_unikmers));
            std::optional<Unikmer> current_min;
            if (_min != std::end(window_unikmers)) {
                current_min = *_min;
            }
            
            if (unikmer_hasher_on_left) {
                unikmer_hasher.reset();
                for (uint16_t i = this->_K - _unikmer_K; i < this->_K; ++i) {
                    unikmer_hasher.eat(*(kmer_window.begin() + i));
                }
                unikmer_hasher_on_left = false;
            }

            const char front = this->kmer_window.front();
            const char ufront = *(kmer_window.begin() + this->_K - _unikmer_K);

            for (auto& symbol : symbols) {
                window_hasher.update(front, symbol);
                unikmer_hasher.update(ufront, symbol);

                Unikmer cur_unikmer(unikmer_hasher.hashvalue);
                if (!current_min || (ukhs_map->query(cur_unikmer) && cur_unikmer < current_min.value())) {
                    hashes.push_back(shift_type(BinnedKmer(window_hasher.hashvalue, cur_unikmer),
                                     symbol));
                } else {
                    hashes.push_back(shift_type(BinnedKmer(window_hasher.hashvalue, current_min.value()),
                                     symbol));
                }
                // update hashers back to root 0
                window_hasher.reverse_update(front, symbol);
                unikmer_hasher.reverse_update(ufront, symbol);
            }

            return hashes;
        }
    };

};

}
}

#endif
