/**
 * (c) Camille Scott, 2019
 * File   : dbg.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */

#ifndef BOINK_DBG_HH
#define BOINK_DBG_HH

#include "boink/meta.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/kmers/kmerclient.hh"
#include "boink/processors.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/storage_types.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/sequences/exceptions.hh"
#include "boink/traversal.hh"

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>


namespace boink {

template <class StorageType,
          class ShifterType>
class dBG : public dBGWalker<dBG<StorageType, ShifterType>> {

public:

    typedef ShifterType                              shifter_type;
    typedef StorageType                              storage_type;

    typedef dBGWalker<dBG<StorageType, ShifterType>> walker_type;
    typedef typename walker_type::extender_type      extender_type;

    typedef hashing::KmerIterator<ShifterType>       kmer_iter_type;
    typedef typename walker_type::alphabet           alphabet;

    typedef typename walker_type::hash_type          hash_type;
    typedef typename walker_type::kmer_type          kmer_type;

    template<bool Dir>
        using shift_type = hashing::ShiftModel<hash_type, Dir>;

    typedef std::pair<std::vector<shift_type<hashing::DIR_LEFT>>,
                      std::vector<shift_type<hashing::DIR_RIGHT>>> shift_pair_type;

    typedef std::pair<std::vector<kmer_type>,
                      std::vector<kmer_type>>             neighbor_pair_type;



protected:

    std::shared_ptr<StorageType> S;

    using walker_type::_K;

public:

    friend walker_type;

    using walker_type::K;

    //dBG(ShifterType& hasher, std::shared_ptr<StorageType> S);
    dBG(std::shared_ptr<StorageType> S, ShifterType& shifter)
        : walker_type(shifter),
          S(S)
    {
    }

    template<typename... Args>
    dBG(std::shared_ptr<StorageType> S, uint16_t K, Args&&... args)
        : walker_type(K, std::forward<Args>(args)...),
          S(S)
    {
    }

    /**
     * @Synopsis Build a dBG instance owned by a shared_ptr. 
     *
     * @Param args  Variadic args to forward on to dBG constructor
     *
     * @Returns     shared_ptr owning the dBG.
     */

    //static std::shared_ptr<dBG<StorageType, ShifterType>>
    //__attribute__((used)) build(ShifterType& hasher, std::shared_ptr<StorageType> S);
    static std::shared_ptr<dBG<StorageType, ShifterType>>
    __attribute__((used)) build(std::shared_ptr<StorageType> S,
                                ShifterType& hasher)  {
        return std::make_shared<dBG<StorageType, ShifterType>>(S, hasher);
    }

    template<typename... Args>
    static std::shared_ptr<dBG<StorageType, ShifterType>>
    __attribute__((used)) build(std::shared_ptr<StorageType> S,
                                uint16_t K,
                                Args&&... args) {
        return std::make_shared<dBG<StorageType, ShifterType>>(S, K, std::forward<Args>(args)...);
    }

    /**
     * @Synopsis  Makes a shallow clone of the dBG.
     *
     * @Returns   shared_ptr owning the clone.
     */
    std::shared_ptr<dBG<StorageType, ShifterType>> clone() {
        return std::make_shared<dBG<StorageType, ShifterType>>(S->clone(),
                                                               *this);
    }

    /**
     * @Synopsis  Hash the k-mer and add it to the dBG storage.
     *
     * @Param kmer The k-mer.
     *
     * @Returns   True if the k-mer was new; false otherwise.
     */
    inline const bool insert(const std::string& kmer) {
        return S->insert(this->hash(kmer).value());
    }

    inline const bool insert(const hash_type& kmer) {
        return S->insert(kmer.value());
    }

    /**
     * @Synopsis  Add the given hash value to the dBG storage and
     *            return its count after inserion.
     *
     * @Param kmer The k-mer.
     *
     * @Returns    Post-insertion count of the element.
     */
    inline const storage::count_t insert_and_query(hash_type& kmer) {
        return S->insert_and_query(kmer.value());
    }

    inline const storage::count_t insert_and_query(const std::string& kmer) {
        return S->insert_and_query(this->hash(kmer).value());
    }

    /**
     * @Synopsis  Gets the count of the given k-mer.
     *
     * @Param kmer The k-mer.
     *
     * @Returns   The count of the k-mer.
     */
    const storage::count_t query(const std::string& kmer) {
        return S->query(this->hash(kmer).value());
    }

    const storage::count_t query(const hash_type& hashed_kmer) const  {
        return S->query(hashed_kmer.value());
    }

    /**
     * @Synopsis  Number of unique k-mers in the storage.
     *
     * @Returns   Number of unique k-mers.
     */
    uint64_t n_unique() const {
        return S->n_unique_kmers();
    }

    /**
     * @Synopsis  Number of occupied buckets in the storage.
     *
     * @Returns   Occupied buckets.
     */
    uint64_t n_occupied() const {
        return S->n_occupied();
    }

    /**
     * @Synopsis  Gets the length K-1 suffix of the given string.
     *
     * @Param kmer The string to take the suffix of.
     *
     * @Returns   The suffix.
     */
    const std::string suffix(const std::string& kmer) {
        return kmer.substr(kmer.length() - this->_K + 1);
    }

    /**
     * @Synopsis  Gets the length K-1 prefix of the given string.
     *
     * @Param kmer The string to take the prefix of.
     *
     * @Returns   The prefix.
     */
    const std::string prefix(const std::string& kmer) {
        return kmer.substr(0, this->_K - 1);
    }

    /**
     * @Synopsis  The address of the raw underling table, if available.
     *
     * @Returns   
     */
    uint8_t ** get_raw() const {
        return S->get_raw_tables();
    }

    /**
     * @Synopsis  If the StorageType is probabilistic, estimate its false positive rate.
     *
     * @Returns   The false-positive rate; 0 if an exact structure.
     */
    template<typename Dummy = double>
    auto estimated_fp() 
    -> std::enable_if_t<storage::is_probabilistic<StorageType>::value, Dummy>
    {
        return S->estimated_fp();
    }

    /**
     * @Synopsis  Insert all k-mers from the given sequence.
     *
     * @Param sequence string containg the sequence, must be length >= K.
     * @Param kmer_hashes Hash values corresponding to the k-mers.
     * @Param counts The k-mer counts after insertion.
     *
     * @Returns Number of new unique k-mers inserted.
     */
    uint64_t insert_sequence(const std::string&          sequence,
                             std::vector<hash_type>&  kmer_hashes,
                             std::vector<storage::count_t>& counts) {
    
        hashing::KmerIterator<ShifterType> iter(sequence, static_cast<ShifterType*>(this));

        uint64_t         n_consumed = 0;
        size_t           pos = 0;
        storage::count_t count;
        while(!iter.done()) {
            auto h = iter.next();
            count  = insert_and_query(h);

            kmer_hashes.push_back(h);
            counts.push_back(count);

            n_consumed += (count == 1);
            ++pos;
        }

        return n_consumed;
    }

    uint64_t insert_sequence(const std::string&      sequence,
                             std::set<hash_type>& new_kmers) {

        hashing::KmerIterator<ShifterType> iter(sequence, this);

        uint64_t n_consumed = 0;
        size_t pos = 0;
        bool is_new;
        while(!iter.done()) {
            hash_type h = iter.next();
            if(insert(h)) {
                new_kmers.insert(h);
                n_consumed += is_new;
            }
            ++pos;
        }

        return n_consumed;
    }

    uint64_t insert_sequence(const std::string& sequence) {
    
        hashing::KmerIterator<ShifterType> iter(sequence, this);

        uint64_t n_consumed = 0;
        while(!iter.done()) {
            hash_type h = iter.next();
            n_consumed += insert(h);
        }

        return n_consumed;
    }

    /**
     * @Synopsis  Insert the sequence and return the post-insertion k-mer counts.
     *
     * @Param sequence String with sequence to insert, >= length K.
     *
     * @Returns  The counts. 
     */
    std::vector<storage::count_t> insert_and_query_sequence(const std::string& sequence)  {

        hashing::KmerIterator<ShifterType> iter(sequence, this);
        std::vector<storage::count_t> counts(sequence.length() - _K + 1);

        size_t pos = 0;
        while(!iter.done()) {
            hash_type h = iter.next();
            counts[pos] = S->insert_and_query(h.value());
            ++pos;
        }

        return counts;
    }

    /**
     * @Synopsis  Queries the k-mer counts for all k-mers in the sequence.
     *
     * @Param sequence String with sequence, >= length K.
     *
     * @Returns   The counts.
     */
    std::vector<storage::count_t> query_sequence(const std::string& sequence)  {

        hashing::KmerIterator<ShifterType> iter(sequence, this);
        std::vector<storage::count_t> counts(sequence.length() - _K + 1);

        size_t pos = 0;
        while(!iter.done()) {
            hash_type h = iter.next();
            counts[pos] = query(h);
            ++pos;
        }

        return counts;
    }

    void query_sequence(const std::string&             sequence,
                        std::vector<storage::count_t>& counts,
                        std::vector<hash_type>&  hashes) {

        hashing::KmerIterator<ShifterType> iter(sequence, this);

        while(!iter.done()) {
            hash_type h = iter.next();
            storage::count_t result = query(h);
            counts.push_back(result);
            hashes.push_back(h);
        }
    }

    void query_sequence(const std::string& sequence,
                        std::vector<storage::count_t>& counts,
                        std::vector<hash_type>& hashes,
                        std::set<hash_type>& new_hashes) {

        hashing::KmerIterator<ShifterType> iter(sequence, this);

        while(!iter.done()) {
            hash_type h = iter.next();
            auto result = query(h);
            if (result == 0) {
                new_hashes.insert(h);
            }
            counts.push_back(result);
            hashes.push_back(h);
        }
    }

    void save(std::string filename) {
        S->save(filename, _K);
    }

    void load(std::string filename) {
        uint16_t ksize = _K;
        S->load(filename, ksize);
    }

    void serialize(std::string& filename) {
        std::ofstream out(filename.c_str(), std::ios::binary);
        out.write(Tagged<ShifterType>::name_string().c_str(),
                  Tagged<ShifterType>::NAME.size());
        out.write(Tagged<ShifterType>::version_binary(),
                  sizeof(Tagged<ShifterType>::OBJECT_ABI_VERSION));

    }

    /**
     * @Synopsis  Remove all k-mers from the dBG / set all values to zero.
     */
    void reset() {
        S->reset();
    }

    auto get_hash_iter(const std::string& sequence)
    -> std::shared_ptr<hashing::KmerIterator<ShifterType>> {

        return std::make_shared<hashing::KmerIterator<ShifterType>>(sequence, this);
    }

    ShifterType get_hasher() {
        return ShifterType(*this);
    }

    using Processor = InserterProcessor<dBG>;

};


extern template class dBG<storage::BitStorage, hashing::FwdRollingShifter>;
extern template class dBG<storage::BitStorage, hashing::CanRollingShifter>;
extern template class dBG<storage::BitStorage, hashing::FwdUnikmerShifter>;
extern template class dBG<storage::BitStorage, hashing::CanUnikmerShifter>;

extern template class dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>;
extern template class dBG<storage::SparseppSetStorage, hashing::CanRollingShifter>;
extern template class dBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>;
extern template class dBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>;

extern template class dBG<storage::ByteStorage, hashing::FwdRollingShifter>;
extern template class dBG<storage::ByteStorage, hashing::CanRollingShifter>;
extern template class dBG<storage::ByteStorage, hashing::FwdUnikmerShifter>;
extern template class dBG<storage::ByteStorage, hashing::CanUnikmerShifter>;

extern template class dBG<storage::NibbleStorage, hashing::FwdRollingShifter>;
extern template class dBG<storage::NibbleStorage, hashing::CanRollingShifter>;
extern template class dBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>;
extern template class dBG<storage::NibbleStorage, hashing::CanUnikmerShifter>;

extern template class dBG<storage::QFStorage, hashing::FwdRollingShifter>;
extern template class dBG<storage::QFStorage, hashing::CanRollingShifter>;
extern template class dBG<storage::QFStorage, hashing::FwdUnikmerShifter>;
extern template class dBG<storage::QFStorage, hashing::CanUnikmerShifter>;

extern template class dBGWalker<dBG<storage::BitStorage, hashing::FwdRollingShifter>>;
extern template class dBGWalker<dBG<storage::BitStorage, hashing::CanRollingShifter>>;
extern template class dBGWalker<dBG<storage::BitStorage, hashing::FwdUnikmerShifter>>;
extern template class dBGWalker<dBG<storage::BitStorage, hashing::CanUnikmerShifter>>;

extern template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>>;
extern template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::CanRollingShifter>>;
extern template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>>;
extern template class dBGWalker<dBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>>;

extern template class dBGWalker<dBG<storage::ByteStorage, hashing::FwdRollingShifter>>;
extern template class dBGWalker<dBG<storage::ByteStorage, hashing::CanRollingShifter>>;
extern template class dBGWalker<dBG<storage::ByteStorage, hashing::FwdUnikmerShifter>>;
extern template class dBGWalker<dBG<storage::ByteStorage, hashing::CanUnikmerShifter>>;

extern template class dBGWalker<dBG<storage::NibbleStorage, hashing::FwdRollingShifter>>;
extern template class dBGWalker<dBG<storage::NibbleStorage, hashing::CanRollingShifter>>;
extern template class dBGWalker<dBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>>;
extern template class dBGWalker<dBG<storage::NibbleStorage, hashing::CanUnikmerShifter>>;

extern template class dBGWalker<dBG<storage::QFStorage, hashing::FwdRollingShifter>>;
extern template class dBGWalker<dBG<storage::QFStorage, hashing::CanRollingShifter>>;
extern template class dBGWalker<dBG<storage::QFStorage, hashing::FwdUnikmerShifter>>;
extern template class dBGWalker<dBG<storage::QFStorage, hashing::CanUnikmerShifter>>;

extern template class hashing::KmerIterator<dBG<storage::BitStorage, hashing::FwdRollingShifter>>;
extern template class hashing::KmerIterator<dBG<storage::BitStorage, hashing::CanRollingShifter>>;
extern template class hashing::KmerIterator<dBG<storage::BitStorage, hashing::FwdUnikmerShifter>>;
extern template class hashing::KmerIterator<dBG<storage::BitStorage, hashing::CanUnikmerShifter>>;

extern template class hashing::KmerIterator<dBG<storage::SparseppSetStorage, hashing::FwdRollingShifter>>;
extern template class hashing::KmerIterator<dBG<storage::SparseppSetStorage, hashing::CanRollingShifter>>;
extern template class hashing::KmerIterator<dBG<storage::SparseppSetStorage, hashing::FwdUnikmerShifter>>;
extern template class hashing::KmerIterator<dBG<storage::SparseppSetStorage, hashing::CanUnikmerShifter>>;

extern template class hashing::KmerIterator<dBG<storage::ByteStorage, hashing::FwdRollingShifter>>;
extern template class hashing::KmerIterator<dBG<storage::ByteStorage, hashing::CanRollingShifter>>;
extern template class hashing::KmerIterator<dBG<storage::ByteStorage, hashing::FwdUnikmerShifter>>;
extern template class hashing::KmerIterator<dBG<storage::ByteStorage, hashing::CanUnikmerShifter>>;

extern template class hashing::KmerIterator<dBG<storage::NibbleStorage, hashing::FwdRollingShifter>>;
extern template class hashing::KmerIterator<dBG<storage::NibbleStorage, hashing::CanRollingShifter>>;
extern template class hashing::KmerIterator<dBG<storage::NibbleStorage, hashing::FwdUnikmerShifter>>;
extern template class hashing::KmerIterator<dBG<storage::NibbleStorage, hashing::CanUnikmerShifter>>;

extern template class hashing::KmerIterator<dBG<storage::QFStorage, hashing::FwdRollingShifter>>;
extern template class hashing::KmerIterator<dBG<storage::QFStorage, hashing::CanRollingShifter>>;
extern template class hashing::KmerIterator<dBG<storage::QFStorage, hashing::FwdUnikmerShifter>>;
extern template class hashing::KmerIterator<dBG<storage::QFStorage, hashing::CanUnikmerShifter>>;


}

#endif
