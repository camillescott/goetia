/**
 * (c) Camille Scott, 2019
 * File   : hashextender.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 08.01.2020
 */

#ifndef BOINK_HASHEXTENDER_HH
#define BOINK_HASHEXTENDER_HH

#include <vector>

#include "boink/hashing/hashshifter.hh"
#include "boink/hashing/canonical.hh"
#include "boink/hashing/kmer_span.hh"
#include "boink/sequences/alphabets.hh"
#include "boink/is_detected.hh"

#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/ukhshashshifter.hh"

namespace boink::hashing {


// Detectors for delegating extension machinery    
template<class ShifterType>
using left_extension_t =
    decltype(std::declval<ShifterType&>().left_extensions());

template<class ShifterType>
using supports_left_extension = is_detected<left_extension_t, ShifterType>;


template<class ShifterType>
using right_extension_t =
    decltype(std::declval<ShifterType&>().right_extensions());

template<class ShifterType>
using supports_right_extension = is_detected<right_extension_t, ShifterType>;




/**
 * @Synopsis  Find left and right extensions from a k-mer. This is stateful,
 *            with the current k-mer stored in a ring buffer. A HashShifter
 *            child class is passed to provide a shifting policy.
 *
 *            Sometimes, the Shifting policy needs to track k-mer state itself;
 *            for example, when computing minimizers along with hashes. The
 *            KmerSpanMixin is templated on the ShifterType, and the mixin
 *            does some template magic to resolve whether ShifterType already
 *            has its own kmer ring buffer. When it does, KmerSpanMixin is
 *            empty, which (sort of hackily) statically resolves the diamond
 *            inheritance problem without virtual inheritance. Given that this
 *            class has some of the hottest code paths in the library, avoiding
 *            vtable lookups is useful.
 *
 * @tparam ShifterType   HashShifter providing the shifting policy.
 */
template <class ShifterType>
class HashExtender : public ShifterType,
                     public KmerSpanMixin<ShifterType>::type {

public:

    typedef typename KmerSpanMixin<ShifterType>::type span_mixin_type;

    typedef ShifterType                      shifter_type;
    typedef typename shifter_type::hash_type hash_type;
    typedef typename hash_type::value_type   value_type;

    template<Direction_t D>
        using shift_type = ShiftModel<hash_type, D>;
    typedef ShiftModel<hash_type, DIR_LEFT>  shift_left_type;
    typedef ShiftModel<hash_type, DIR_RIGHT> shift_right_type;

    typedef KmerModel<hash_type>             kmer_type;

    using ShifterType::_K;

    typedef typename ShifterType::alphabet alphabet;

    template<typename... ExtraArgs>
    HashExtender(const std::string& start,
                 uint16_t           K,
                 ExtraArgs&&...     args)
        : ShifterType(K, std::forward<ExtraArgs>(args)...),
          span_mixin_type(K)
    {
        set_cursor(start);
    }

    template<typename... ExtraArgs>
    HashExtender(uint16_t K,
                 ExtraArgs&&... args)
        : ShifterType(K, std::forward<ExtraArgs>(args)...),
          span_mixin_type(K)
    {
    }

    HashExtender(const HashExtender& extender)
        : ShifterType(static_cast<ShifterType>(extender)),
          span_mixin_type(extender.K())
    {
    }

    HashExtender(const shifter_type& shifter)
        : ShifterType(shifter),
          span_mixin_type(shifter.K())
    {
    }

    using shifter_type::get;
    using shifter_type::shift_right;
    using shifter_type::shift_left;
    using shifter_type::K;


    /**
     * @Synopsis  Shift cursor left from current value using
     *            symbol c and return hash value.
     *
     * @Param c
     *
     * @Returns   
     */
    hash_type shift_left(const char& c) {
        if (!this->is_loaded()) {
            throw UninitializedShifterException();
        }

        hash_type h = shift_left(c, this->back());
        this->kmer_window.push_front(c);
        return h;
    }

    /**
     * @Synopsis  Gather the left extensions from the current position using
     *            the given symbols, the default being those from our alphabet.
     *
     *            NOTE: there is template wonk here. This version is called by
     *            default, that is, ShifterType is a normal CRTP derivation
     *            from HashShifter and does not have a _left_extensions
     *            member function.
     *
     * @tparam Dummy    SFINAE dummy var. Do not specialize!
     * @Param symbols   Alphabet to extend with.
     *
     * @Returns   A list of left extensions.
     */
    template<typename Dummy = std::vector<shift_left_type>>
    auto left_extensions(const std::string& symbols = alphabet::SYMBOLS)
    -> std::enable_if_t<!supports_left_extension<ShifterType>::value, Dummy> {

        if (!this->is_loaded()) {
            throw UninitializedShifterException();
        }

        std::vector<shift_left_type> hashes;
        auto back = this->back();
        for (auto symbol : symbols) {
            hash_type h = this->shift_left(symbol, back);
            shift_left_type result(h, symbol);
            hashes.push_back(result);
            this->shift_right(symbol, back);
        }

        return hashes;
    }

    /**
     * @Synopsis  Gather the left extensions from the current position using
     *            the given symbols, the default being those from our alphabet.
     *
     *            NOTE: template wonk! If ShifterType has its own _left_extensions
     *            member function, we delegate to that instead. This is for when
     *            it's slow to shift left-right-...-left-right repeatedly and
     *            the shifter needs to provide its own implementation. This also
     *            avoids a ton of duplicated code via partial template specialization.
     *
     * @tparam Dummy    SFINAE dummy var. Do not specialize!
     * @Param symbols   Alphabet to extend with.
     *
     * @Returns   A list of left extensions.
     */
    template<typename Dummy = std::vector<shift_left_type>>
    auto left_extensions(const std::string& symbols = alphabet::SYMBOLS)
    -> std::enable_if_t<supports_left_extension<ShifterType>::value, Dummy> {

        if (!this->is_loaded()) {
            throw UninitializedShifterException();
        }

        return ShifterType::_left_extensions(symbols);
    }

    /**
     * @Synopsis  Shift cursor right from current value using
     *            symbol c and return hash value.
     *
     * @Param c
     *
     * @Returns   
     */
    hash_type shift_right(const char& c) {
        if (!this->is_loaded()) {
            throw UninitializedShifterException();
        }

        hash_type h = shift_right(this->front(), c);
        this->kmer_window.push_back(c);
        return h;
    }

    /**
     * @Synopsis  Gather the right extensions from the current position using
     *            the given symbols, the default being those from our alphabet.
     *
     *            NOTE: there is template wonk here. This version is called by
     *            default, that is, ShifterType is a normal CRTP derivation
     *            from HashShifter and does not have a _left_extensions
     *            member function.
     *
     * @tparam Dummy    SFINAE dummy var. Do not specialize!
     * @Param symbols   Alphabet to extend with.
     *
     * @Returns   A vector of right extensions.
     */
    template<typename Dummy = std::vector<shift_right_type>>
    auto right_extensions(const std::string& symbols = alphabet::SYMBOLS)
    -> std::enable_if_t<!supports_right_extension<ShifterType>::value, Dummy> {

        if (!this->is_loaded()) {
            throw UninitializedShifterException();
        }

        std::vector<shift_right_type> hashes;
        auto front = this->front();
        for (auto symbol : symbols) {
            hash_type h = this->shift_right(front, symbol);
            hashes.push_back(shift_right_type(h, symbol));
            this->shift_left(front, symbol);
        }
        return hashes;
    }

    /**
     * @Synopsis  Gather the right extensions from the current position using
     *            the given symbols, the default being those from our alphabet.
     *
     *            NOTE: template wonk! If ShifterType has its own _left_extensions
     *            member function, we delegate to that instead. This is for when
     *            it's slow to shift left-right-...-left-right repeatedly and
     *            the shifter needs to provide its own implementation. This also
     *            avoids a ton of duplicated code via partial template specialization.
     *
     * @tparam Dummy    SFINAE dummy var. Do not specialize!
     * @Param symbols   Alphabet to extend with.
     *
     * @Returns   A vector of right extensions.
     */
    template<typename Dummy = std::vector<shift_right_type>>
    auto right_extensions(const std::string& symbols = alphabet::SYMBOLS)
    -> std::enable_if_t<supports_right_extension<ShifterType>::value, Dummy> {

        if (!this->is_loaded()) {
            throw UninitializedShifterException();
        }

        return ShifterType::_right_extensions(symbols);
    }

    /**
     * @Synopsis  Initialize the base hash function and set the kmer span cursor.
     *
     * @Param sequence
     *
     * @Returns   
     */
    hash_type set_cursor(const std::string& sequence) {
        this->hash_base(sequence.c_str());
        this->load(sequence);
        return get();
    }

    /**
     * @Synopsis  Initialize the base hash function and set the kmer span cursor.
     *
     * @Param sequence
     *
     * @Returns   
     */
    hash_type set_cursor(const char * sequence) {
        // less safe! does not check length
        this->hash_base(sequence);
        this->load(sequence);
        return get();
    }

    /**
     * @Synopsis  Get the cursor as a string.
     *
     * @Returns   
     */
    const std::string get_cursor() const {
        return this->to_string();
    }

    /**
     * @Synopsis  Load the cursor on to a deque<char>
     *
     * @Param d
     */
    void get_cursor(std::deque<char>& d) const {
        this->to_deque(d);
    }

    /**
     * @Synopsis  Get the current cursor and hash value as a kmer_type.
     *
     * @Param result
     */
    void get_cursor(kmer_type& result) {
        result.hash = get();
        result.kmer = get_cursor();
    }
};


}

#endif
