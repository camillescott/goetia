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
#include <limits>

#include "boink/boink.hh"
#include "boink/meta.hh"
#include "boink/sequences/alphabets.hh"
#include "boink/storage/sparsepp/spp.h"

#include "boink/hashing/canonical.hh"
#include "boink/hashing/kmer_span.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/hashshifter.hh"


namespace boink::hashing { 

void test();

template <class ShifterType>
struct UKHS {

    typedef typename ShifterType::hash_type hash_type;
    typedef typename hash_type::value_type  value_type;

    typedef Partitioned<hash_type>          Unikmer;


protected:

    // map of unikmer hashes to partitions
    typedef spp::sparse_hash_map<value_type, uint64_t> pmap_t;
    pmap_t                                             pmap;
    std::vector<hash_type>                             hashes;


public:

    // Window size W; will correspond to K in other dBG contexts
    const uint16_t W;
    // The minimizer length
    const uint16_t K;

    explicit UKHS(uint16_t W,
                 uint16_t K,
                 std::vector<std::string>& ukhs) 
        : K (K),
          W (W)
    {
        if (ukhs.front().size() != K) {
            throw BoinkException("K does not match k-mer size from provided UKHS");
        }

        uint64_t pid = 0;
        for (const auto& unikmer : ukhs) {
            hash_type h = ShifterType::hash(unikmer, K);

            if (!pmap.count(h.value())) {
                hashes.push_back(h);
                pmap[h.value()] = pid;
                ++pid;
            }
        }
    }

    static std::shared_ptr<UKHS> build(uint16_t W,
                                      uint16_t K,
                                      std::vector<std::string>& ukhs) {
        return std::make_shared<UKHS>(W, K, ukhs);
    }

    std::optional<Unikmer> query(hash_type unikmer_hash) {

        auto search = pmap.find(unikmer_hash.value());
        if (search != pmap.end()) {
            return {{unikmer_hash, search->second}};
        } else {
            return {};
        }
    }

    std::vector<hash_type> get_hashes() const {
        return hashes;
    }

    const size_t n_hashes() const {
        return pmap.size();
    }
};

template<typename T>
struct UnikmerShifterPolicy;

template<template <typename, typename> typename ShiftPolicy,
                                       typename HashType,
                                       typename Alphabet>
class UnikmerShifterPolicy<ShiftPolicy<HashType, Alphabet>> : public KmerSpan {

    /* Shifter that keeps track of each k-mer's associated
     * Unikmer. The hash_type will end being composed as:
     *
     *     hash_type:
     *         uint64_t              hash
     *         Partitioned<uint64_t> minimizer:
     *             uint64_t value
     *             size_t   partition
     * 
     * At least, assuming value_type is uint64_t, which it basically
     * always is. That is, hash_type follows WmerModel<uint64_t, Unikmer>, 
     * and kmer_type will have WmerModel as its hash.
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
     * 
     * NOTE: This is both a valid shifter policy for HashShifter AND
     *       and valid extension policy for HashExtender.
     */


public:

    // policy typedefs for HashShifter
    typedef ShiftPolicy<HashType, Alphabet>                     base_shift_policy;
    typedef HashShifter<base_shift_policy>                      base_shifter_type;

    typedef UnikmerShifterPolicy<base_shift_policy>             type;

    // Shared typedefs
    typedef UKHS<base_shifter_type>                             ukhs_type;
    typedef typename ukhs_type::Unikmer                         minimizer_type;
    typedef WmerModel<typename ukhs_type::hash_type,
                      typename ukhs_type::Unikmer>              hash_type;
    typedef KmerModel<hash_type>                                kmer_type;
    typedef typename hash_type::value_type                      value_type;
    typedef Alphabet                                            alphabet;

    typedef hash_type                                           wmer_type;
    typedef ShiftModel<wmer_type, DIR_LEFT>                     shift_left_type;
    typedef ShiftModel<wmer_type, DIR_RIGHT>                    shift_right_type;

    static constexpr bool has_kmer_span = true;

protected:

    // K size of the unikmer
    uint16_t _unikmer_K;

    base_shifter_type window_hasher;
    base_shifter_type unikmer_hasher;

    // Current position of the unikmer hasher within
    // the window: it will need to be moved across the window
    // if we change directions
    bool unikmer_hasher_on_left;

    // We'll store the unikmers for the current window and their
    // indices, only finding the minimum when queried.
    std::deque<minimizer_type> window_unikmers;
    std::deque<size_t>         unikmer_indices;

    void shift_unikmers_left() {
        for (auto& index : unikmer_indices) {
            ++index;
        }
        if (unikmer_indices.back() > K - _unikmer_K) {
            unikmer_indices.pop_back();
            window_unikmers.pop_back();
        }
    }

    void update_unikmer_left(const char c) {
        if (!unikmer_hasher_on_left) {
            // if we're not on the left of the window, reset the cursor there
            unikmer_hasher.hash_base(ring.begin(),
                                     ring.begin() + _unikmer_K);

            unikmer_hasher_on_left = true;
        }
        // othewise just shift the new symbol on
        // std::cout << "Reverse update, "
        //          << *(ring.begin() + _unikmer_K - 1)
        //          << " => "
        //          << c << std::endl;
        unikmer_hasher.shift_left(c, *(ring.begin() + _unikmer_K - 1));

        shift_unikmers_left();
        auto unikmer = ukhs_map->query(unikmer_hasher.get());
        if (unikmer) {
            unikmer_indices.push_front(0);
            window_unikmers.push_front(unikmer.value());
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
        //std::cout << "UnikmerShifter update_unikmer_right(" 
        //          << c << "): " << this->get_cursor() << std::endl;
        // if the ukhs hasher is on the left, reset the cursor
        if (unikmer_hasher_on_left) {
            //std::cout << "UnikmerShifter set_unikmer_hasher_right_suffix" << std::endl;
            unikmer_hasher.hash_base(ring.begin() + K - _unikmer_K,
                                     ring.end());
            unikmer_hasher_on_left = false;
        }
            //std::cout << "UnikmerShifter update: "
            //          << *(ring.begin() + K - _unikmer_K)
            //          << " => " << c << std::endl;
        unikmer_hasher.shift_right(*(ring.begin() + K - _unikmer_K), c);
        
        shift_unikmers_right();
        auto unikmer = ukhs_map->query(unikmer_hasher.get());
        if (unikmer) {
            //std::cout << "UnikmerShifter: found unikmer" << std::endl;
            unikmer_indices.push_back(K - _unikmer_K);
            window_unikmers.push_back(unikmer.value());
        }
        //std::cout << "UnikmerShifter: leave update_unikmer, uhash: " << cur_unikmer << std::endl;
    }

    static minimizer_type get_min_unikmer(std::deque<minimizer_type>& window_unikmers) {
        //std::cout << "UnikmerShifter: " << repr(window_unikmers) << std::endl;
        //std::cout << "UnikmerShifter: " << repr(unikmer_indices) << std::endl;

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

    // W,K unikmer map
    std::shared_ptr<ukhs_type> ukhs_map;
    const uint16_t K;

    uint16_t unikmer_K() const {
        return _unikmer_K;
    }

    wmer_type hash_base_impl(const char * sequence) {
        wmer_type h = _hash(sequence,
                            window_hasher,
                            unikmer_hasher,
                            ukhs_map,
                            window_unikmers,
                            unikmer_indices);
        this->load(sequence);
        unikmer_hasher_on_left = false;

        return h;
        //std::cout << "UnikmerShifter: finish _hash_base(): " << this->get() << std::endl;
    }

    template<class It>
    wmer_type hash_base_impl(It begin, It end) {
        // TODO: this does extra copying
        this->load(begin, end);
        wmer_type h = _hash(this->to_string().c_str(),
                            window_hasher,
                            unikmer_hasher,
                            ukhs_map,
                            window_unikmers,
                            unikmer_indices);
        unikmer_hasher_on_left = false;

        return h;
    }

    wmer_type get_impl() {
        typename base_shifter_type::hash_type hash = window_hasher.get();
        minimizer_type minimizer = get_min_unikmer(window_unikmers);
        return {hash, minimizer};
    }

    static wmer_type _hash(const char * sequence, 
                           const uint16_t W,
                           const uint16_t K,
                           const std::shared_ptr<ukhs_type>& ukhs) {

        UnikmerShifterPolicy shifter(W,
                               K,
                               ukhs);
        return shifter._hash(sequence,
                             shifter.window_hasher,
                             shifter.unikmer_hasher,
                             shifter.ukhs_map,
                             shifter.window_unikmers,
                             shifter.unikmer_indices);
    }

    static wmer_type _hash(const char *                sequence,
                           base_shifter_type&               window_hasher,
                           base_shifter_type&               unikmer_hasher,
                           std::shared_ptr<ukhs_type>& ukhs_map,
                           std::deque<minimizer_type>& window_unikmers,
                           std::deque<size_t>&         unikmer_indices) {

        window_hasher.hash_base(sequence);
        unikmer_hasher.hash_base(sequence);
        
        unikmer_indices.clear();
        window_unikmers.clear();

        auto unikmer = ukhs_map->query(unikmer_hasher.get());
        if (unikmer) {
            unikmer_indices.push_back(0);
            window_unikmers.push_back(unikmer.value());
        }

        for (uint16_t i = unikmer_hasher.K; i < window_hasher.K; ++i) {
            // once we've eaten the first K bases, start keeping
            // track of unikmers

            unikmer_hasher.shift_right(sequence[i - unikmer_hasher.K],
                                       sequence[i]);

            unikmer = ukhs_map->query(unikmer_hasher.get());
            if (unikmer) {
                unikmer_indices.push_back(i - unikmer_hasher.K + 1);
                window_unikmers.push_back(unikmer.value());
            }

        }

        return {window_hasher.get(), get_min_unikmer(window_unikmers)};
    }


    wmer_type shift_left_impl(const char in, const char out) {
        // update main window left
        window_hasher.shift_left(in, out);
        update_unikmer_left(in);
        ring.push_front(in);
        return get_impl();
    }

    wmer_type shift_left_impl(const char in) {
        // update main window left
        window_hasher.shift_left(in, this->back());
        update_unikmer_left(in);
        ring.push_front(in);
        return get_impl();
    }

    wmer_type shift_right_impl(const char out, const char in) {
        // update main window right
        window_hasher.shift_right(out, in);
        update_unikmer_right(in);
        ring.push_back(in);
        return get_impl();
    }

    wmer_type shift_right_impl(const char in) {
        // update main window right
        window_hasher.shift_right(this->front(), in);
        update_unikmer_right(in);
        ring.push_back(in);
        return get_impl();
    }

    auto left_extensions_impl()
    -> std::vector<shift_left_type> {

        std::vector<shift_left_type> hashes;

        // First get the min unikmer in the W-1 prefix, if there is one
        auto _last = unikmer_indices.back() > K - _unikmer_K ? 
                     std::end(window_unikmers) - 1 :
                     std::end(window_unikmers);
        auto _min = std::min_element(std::begin(window_unikmers),
                                     _last);
        std::optional<minimizer_type> current_min;
        if (_min != _last) {
            current_min = *_min;
        }
        
        if (!unikmer_hasher_on_left) {
            unikmer_hasher.hash_base(ring.begin(),
                                     ring.begin() + _unikmer_K);

            unikmer_hasher_on_left = true;
        }

        // now compare each neighbor unikmer to the current min
        // to see if the neighbor w-mer has a different minimizer
        const char back = this->ring.back();
        const char uback = *(this->ring.begin() + _unikmer_K - 1);
        for (const auto& symbol : alphabet::SYMBOLS) {
            window_hasher.shift_left(symbol, back);
            unikmer_hasher.shift_left(symbol, uback);

            auto unikmer = ukhs_map->query(unikmer_hasher.get());
            if (!current_min || (unikmer && unikmer.value() < current_min.value())) {
                hashes.emplace_back(shift_left_type{wmer_type{window_hasher.get(), unikmer.value()},
                                                    symbol});
            } else {
                hashes.emplace_back(shift_left_type{wmer_type{window_hasher.get(), current_min.value()},
                                                    symbol});
            }

            // update hashers back to root 0
            window_hasher.shift_right(symbol, back);
            unikmer_hasher.shift_right(symbol, uback);

        }

        return hashes;
    }

    auto right_extensions_impl()
    -> std::vector<shift_right_type> {
        
        std::vector<shift_right_type> hashes;

        // First get the min unikmer in the W-1 prefix, if there is one
        auto _first = unikmer_indices.front() == 0 ? 
                      std::begin(window_unikmers) + 1 :
                      std::begin(window_unikmers);
        auto _min = std::min_element(_first, std::end(window_unikmers));
        std::optional<minimizer_type> current_min;
        if (_min != std::end(window_unikmers)) {
            current_min = *_min;
        }
        
        if (unikmer_hasher_on_left) {
            unikmer_hasher.hash_base(ring.begin() + K - _unikmer_K,
                                     ring.end());
            unikmer_hasher_on_left = false;
        }

        const char front = this->ring.front();
        const char ufront = *(ring.begin() + K - _unikmer_K);

        for (const auto& symbol : alphabet::SYMBOLS) {
            window_hasher.shift_right(front, symbol);
            unikmer_hasher.shift_right(ufront, symbol);

            auto unikmer = ukhs_map->query(unikmer_hasher.get());
            if (!current_min || (unikmer && unikmer.value() < current_min.value())) {
                hashes.emplace_back(shift_right_type{wmer_type{window_hasher.get(), unikmer.value()},
                                                     symbol});
            } else {
                hashes.emplace_back(shift_right_type{wmer_type{window_hasher.get(), current_min.value()},
                                                     symbol});
            }
            window_hasher.shift_left(front, symbol);
            unikmer_hasher.shift_left(ufront, symbol);
        }

        return hashes;
    }

    hash_type set_cursor_impl(const char * sequence) {
        this->load(sequence);
        hash_base_impl(sequence);
        return get_impl();
    }

protected:

    explicit UnikmerShifterPolicy(uint16_t K,
                                  uint16_t unikmer_K,
                                  std::shared_ptr<ukhs_type> ukhs)
        :  KmerSpan       (K),
           window_hasher  (K),
           unikmer_hasher (unikmer_K),
           K              (K),
           _unikmer_K     (unikmer_K),
           ukhs_map       (std::move(ukhs))
    {
        if (ukhs_map->W != K) {
            throw BoinkException("Shifter K does not match UKHS::Map W.");
        }
        if (ukhs_map->K != unikmer_K) {
            throw BoinkException("Shifter unikmer_K does not match UKHS::Map K.");
        }
    }

    explicit UnikmerShifterPolicy(UnikmerShifterPolicy& other)
        : UnikmerShifterPolicy(other.K,
                               other.unikmer_K(),
                               other.ukhs_map)
    {
    }

    explicit UnikmerShifterPolicy(const UnikmerShifterPolicy& other)
        : UnikmerShifterPolicy(other.K,
                         other.unikmer_K(),
                         other.ukhs_map)
    {
    }
};


/**
 * @brief Works as an alias to convert UnikmerShifter to the HashShifter's
 *        shifting policy template:w
 *  specification.
 * 
 * @tparam HashType 
 * @tparam Alphabet 
 */
template <typename HashType, typename Alphabet=DNA_SIMPLE>
struct UnikmerLemirePolicy : public UnikmerShifterPolicy<RollingHashShifter<HashType, Alphabet>>
{
    protected:
        using UnikmerShifterPolicy<RollingHashShifter<HashType, Alphabet>>::UnikmerShifterPolicy;
};

typedef UnikmerLemirePolicy<HashModel<uint64_t>> FwdUnikmerPolicy;
typedef UnikmerLemirePolicy<CanonicalModel<uint64_t>> CanUnikmerPolicy;

extern template class UKHS<FwdRollingShifter>;
extern template class UKHS<CanRollingShifter>;

extern template class UnikmerShifterPolicy<FwdLemirePolicy>;
extern template class UnikmerShifterPolicy<CanLemirePolicy>;

extern template class UnikmerLemirePolicy<HashModel<uint64_t>, DNA_SIMPLE>;
extern template class UnikmerLemirePolicy<CanonicalModel<uint64_t>, DNA_SIMPLE>;

}

#endif
