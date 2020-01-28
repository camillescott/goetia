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
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/shifter_types.hh"
#include "boink/sequences/alphabets.hh"

#include "boink/is_detected.hh"



namespace boink::hashing {


// Detectors for delegating extension machinery    
template<class ShifterType>
using left_extension_t =
    decltype(std::declval<ShifterType&>().left_extensions_impl());

template<class ShifterType>
using supports_left_extension = is_detected<left_extension_t, ShifterType>;


template<class ShifterType>
using right_extension_t =
    decltype(std::declval<ShifterType&>().right_extensions_impl());

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

;

template <typename ExtensionPolicy>
class HashExtender : public ExtensionPolicy {

public:

    typedef ExtensionPolicy                         extension_policy;
    typedef typename extension_policy::shifter_type shifter_type;
    typedef typename extension_policy::value_type   value_type;
    typedef typename extension_policy::hash_type    hash_type;
    typedef typename extension_policy::kmer_type    kmer_type;
    typedef typename extension_policy::alphabet     alphabet;

    typedef ShiftModel<hash_type, DIR_LEFT>         shift_left_type;
    typedef ShiftModel<hash_type, DIR_RIGHT>        shift_right_type;

    using extension_policy::K;

    template<typename... ExtraArgs>
    explicit HashExtender(const std::string& start,
                          uint16_t           K,
                          ExtraArgs&&...     args)
        : extension_policy(start, K, std::forward<ExtraArgs>(args)...)
    {
        //std::cout << "END HashExtender(start, K...) ctor " << this << " / " << static_cast<extension_policy*>(this) << std::endl;
    }

    template<typename... ExtraArgs>
    explicit HashExtender(uint16_t K,
                          ExtraArgs&&... args)
        : extension_policy(K, std::forward<ExtraArgs>(args)...)
    {
        //std::cout << "END HashExtender(K...) ctor " << this << " / " << static_cast<extension_policy*>(this) << std::endl;
    }

    explicit HashExtender(const HashExtender& extender)
        : extension_policy(static_cast<const extension_policy &>(extender))
    {
        //std::cout << "END HashExtender(HashExtender...) ctor " << this << " / " << static_cast<extension_policy*>(this) << std::endl;
    }

    explicit HashExtender(const shifter_type& shifter)
        : extension_policy(shifter)
    {
        //std::cout << "END HashExtender(shifter_type...) ctor " << this << " / " << static_cast<extension_policy*>(this) << std::endl;
    }

    HashExtender() = delete;

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

        return this->shift_left_impl(c);
    }

    hash_type shift_left(const char& in, const char& out) {
        return shift_left(in);
    }

    /**
     * @Synopsis  Gather the left extensions from the current position using
     *            the given symbols, the default being those from our alphabet.
     *
     * @Param symbols   Alphabet to extend with.
     *
     * @Returns   A list of left extensions.
     */
    std::vector<shift_left_type> left_extensions() {

        if (!this->is_loaded()) {
            throw UninitializedShifterException();
        }

        return this->left_extensions_impl();
    }

    /**
     * @Synopsis  Shift cursor right from current value using
     *            symbol c and return hash value.
     *
     * @Param c
     *
     * @Returns   
     */
    hash_type shift_right(const char& c){
        if (!this->is_loaded()) {
            throw UninitializedShifterException();
        }

        return this->shift_right_impl(c);

    }

    hash_type shift_right(const char& out, const char& in) {
        return shift_right(in);
    }

    /**
     * @Synopsis  Gather the right extensions from the current position using
     *            the given symbols, the default being those from our alphabet.
     *
     * @Param symbols   Alphabet to extend with.
     *
     * @Returns   A vector of right extensions.
     */
    std::vector<shift_right_type> right_extensions() {

        if (!this->is_loaded()) {
            throw UninitializedShifterException();
        }

        return this->right_extensions_impl();
    }

    /**
     * @Synopsis  Initialize the base hash function and set the kmer span cursor.
     *
     * @Param sequence
     *
     * @Returns   
     */
    hash_type set_cursor(const std::string& sequence) {
        if (sequence.length() < K) {
            throw SequenceLengthException("Sequence must at least length K");
        }
        return set_cursor(sequence.c_str());
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
        return this->set_cursor_impl(sequence);
    }

    hash_type hash_base(const std::string& sequence) {
        if (sequence.length() < K) {
            throw SequenceLengthException("Sequence must at least length K");
        }

        return set_cursor(sequence);
    }

    hash_type hash_base(const char * sequence) {
        return set_cursor(sequence);
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
        result.hash = this->get();
        result.kmer = get_cursor();
    }
};

template <class ShifterType>
class DefaultExtensionPolicy : public ShifterType,
                               public KmerSpan {

public:

    typedef ShifterType                      shifter_type;
    typedef typename shifter_type::hash_type hash_type;
    typedef typename hash_type::value_type   value_type;

    template<bool Dir>
        using shift_type = ShiftModel<hash_type, Dir>;

    typedef ShiftModel<hash_type, DIR_LEFT>  shift_left_type;
    typedef ShiftModel<hash_type, DIR_RIGHT> shift_right_type;
    typedef KmerModel<hash_type>             kmer_type;
    typedef typename ShifterType::alphabet   alphabet;

    using shifter_type::K;
    using shifter_type::hash;
    using shifter_type::get;
    using shifter_type::is_initialized;

private:

    using shifter_type::shift_right;
    using shifter_type::shift_left;
    using shifter_type::hash_base;

public:


    template<typename... ExtraArgs>
    explicit DefaultExtensionPolicy(uint16_t       K,
                                    ExtraArgs&&... args)
        : ShifterType(K, std::forward<ExtraArgs>(args)...),
          KmerSpan(K)
    {
        //std::cout << "END DefaultExtPolicy(K...) ctor " << this << " / " << static_cast<ShifterType*>(this) << std::endl;
    }

    template<typename... ExtraArgs>
    explicit DefaultExtensionPolicy(const std::string& start,
                                    uint16_t           K,
                                    ExtraArgs&&...     args)
        : ShifterType(start, K, std::forward<ExtraArgs>(args)...),
          KmerSpan(K)
    {
        this->load(start);
        //std::cout << "END DefaultExtPolicy (string&, ...) ctor " << this << " / " << static_cast<ShifterType*>(this) << std::endl;

    }

    explicit DefaultExtensionPolicy(const DefaultExtensionPolicy& policy)
        : ShifterType(static_cast<const ShifterType &>(policy)),
          KmerSpan(policy.K)
    {
        //std::cout << "END DefaultExtPolicy(HashExtender&......) ctor" << this << " / " << static_cast<ShifterType*>(this) << std::endl;
    }

    explicit DefaultExtensionPolicy(const shifter_type& shifter)
        : ShifterType(shifter),
          KmerSpan(shifter.K)
    {
        //std::cout << "END DefaultExtPolicy(shifter_type&...) ctor"   << this << " / " << static_cast<ShifterType*>(this) << std::endl;
    }

    DefaultExtensionPolicy() = delete;

    /**
     * @Synopsis  Shift cursor left from current value using
     *            symbol c and return hash value.
     *
     * @Param c
     *
     * @Returns   
     */
    hash_type shift_left_impl(const char& c) {
        auto h = ShifterType::shift_left(c, this->back());
        this->ring.push_front(c);
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
    std::vector<shift_left_type> left_extensions_impl() {

        std::vector<shift_left_type> hashes;
        auto back = this->back();
        for (const auto& symbol : alphabet::SYMBOLS) {
            hash_type h = ShifterType::shift_left(symbol, back);
            shift_left_type result(h, symbol);
            hashes.push_back(result);
            ShifterType::shift_right(symbol, back);
        }

        return hashes;
    }

    /**
     * @Synopsis  Shift cursor right from current value using
     *            symbol c and return hash value.
     *
     * @Param c
     *
     * @Returns   
     */
    hash_type shift_right_impl(const char& c) {
        auto h = ShifterType::shift_right(this->front(), c);
        this->ring.push_back(c);
        return h;
    }

    /**
     * @Synopsis  Gather the right extensions from the current position using
     *            the given symbols, the default being those from our alphabet.
     *
     * @Param symbols   Alphabet to extend with.
     *
     * @Returns   A vector of right extensions.
     */
     std::vector<shift_right_type> right_extensions_impl() {

        std::vector<shift_right_type> hashes;
        auto front = this->front();
        for (const auto& symbol : alphabet::SYMBOLS) {
            hash_type h = ShifterType::shift_right(front, symbol);
            hashes.push_back(shift_right_type(h, symbol));
            ShifterType::shift_left(front, symbol);
        }
        return hashes;
    }

    hash_type set_cursor_impl(const char * sequence) {
        auto h = ShifterType::hash_base(sequence);
        this->load(sequence);
        return h;
    }
};

typedef HashExtender<DefaultExtensionPolicy<FwdRollingShifter>> FwdRollingExtender;
typedef HashExtender<DefaultExtensionPolicy<CanRollingShifter>> CanRollingExtender;

typedef HashExtender<FwdUnikmerShifter> FwdUnikmerExtender;
typedef HashExtender<CanUnikmerShifter> CanUnikmerExtender;


template <typename ShifterType>
struct extender_selector {
    typedef HashExtender<DefaultExtensionPolicy<ShifterType>> type;
    typedef ShifterType                                       shifter_type;
};


template<>
struct extender_selector<FwdUnikmerShifter> {
    typedef HashExtender<FwdUnikmerShifter> type;
    typedef FwdUnikmerShifter               shifter_type;
};


template<>
struct extender_selector<CanUnikmerShifter> {
    typedef HashExtender<CanUnikmerShifter> type;
    typedef CanUnikmerShifter               shifter_type;
};


template<typename ShifterType>
    using extender_selector_t = typename extender_selector<ShifterType>::type;


extern template class DefaultExtensionPolicy<FwdRollingShifter>;
extern template class DefaultExtensionPolicy<CanRollingShifter>;

extern template class HashExtender<DefaultExtensionPolicy<FwdRollingShifter>>;
extern template class HashExtender<DefaultExtensionPolicy<CanRollingShifter>>;

extern template class HashExtender<FwdUnikmerShifter>;
extern template class HashExtender<CanUnikmerShifter>;

extern template class KmerIterator<FwdRollingExtender>;
extern template class KmerIterator<CanRollingExtender>;

extern template class KmerIterator<FwdUnikmerExtender>;
extern template class KmerIterator<CanUnikmerExtender>;



}

#endif
