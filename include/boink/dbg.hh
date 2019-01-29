/* dbg.hh -- Generic de Bruijn graph
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_DBG_HH
#define BOINK_DBG_HH


#include "boink/assembly.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/storage/storage.hh"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace boink {


template <class StorageType,
          class HashShifter>
class dBG : public hashing::KmerClient,
            public std::enable_shared_from_this<dBG<StorageType, HashShifter>> {

    StorageType S;
    HashShifter hasher;

public:

    const uint16_t N;
    const uint64_t max_table;
    std::vector<uint64_t> sizes;

    typedef HashShifter shifter_type;
	typedef AssemblerMixin<dBG<StorageType, HashShifter>> assembler_type;
    typedef hashing::KmerIterator<HashShifter> kmer_iter_type;
    

    explicit dBG(uint16_t K, const std::vector<uint64_t>& storage_size)
        : KmerClient(K),
          S(storage_size),
          hasher(K),
          sizes(storage_size),
          N(storage_size.size()),
          max_table(*std::max_element(storage_size.begin(), storage_size.end()))
    {    
    }

    explicit dBG(uint16_t K, uint64_t max_table, uint16_t N)
        : KmerClient(K),
          hasher(K),
          sizes(storage::get_n_primes_near_x(N, max_table)),
          N(N),
          max_table(max_table),
          S(storage::get_n_primes_near_x(N, max_table)) 
    {
    }

    std::shared_ptr<dBG<StorageType, HashShifter>> clone() const {
        return std::make_shared<dBG<StorageType, HashShifter>>(_K, sizes);
    }

    /**
     * @Synopsis  Hash a k-mer using the templated HashShifter.
     *
     * @Param kmer String of length K.
     *
     * @Returns   The hash value.
     */
    hashing::hash_t hash(const std::string& kmer) const {
        return hasher.hash(kmer);
    }

    /**
     * @Synopsis  Hash a k-mer using the templated HashShifter.
     *
     * @Param kmer c-string with the k-mer; if longer than K, will only use
     *             first K characters.
     *
     * @Returns   The hash value.
     */
    hashing::hash_t hash(const char * kmer) const {
        return hasher.hash(kmer);
    }

    /**
     * @Synopsis  Hash the k-mer and add it to the dBG storage.
     *
     * @Param kmer The k-mer string.
     *
     * @Returns   True if the k-mer was new; fals otherwise.
     */
    bool add(const std::string& kmer) {
        return S.add(hash(kmer));
    }

    /**
     * @Synopsis  Add the given hash value to the dBG storage.
     *
     * @Param kmer The hash value.
     *
     * @Returns   True if the hash value was new; false otherwise.
     */
    bool add(hashing::hash_t kmer) {
        return S.add(kmer);
    }

    /**
     * @Synopsis  Gets the count of the given k-mer.
     *
     * @Param kmer String with the k-mer.
     *
     * @Returns   The count of the k-mer.
     */
    const storage::count_t get(const std::string& kmer) const {
        return S.get_count(hash(kmer));
    }

    /**
     * @Synopsis  Gets the count of the given hash value.
     *
     * @Param hashed_kmer The hash value.
     *
     * @Returns   The count of the hash value.
     */
    const storage::count_t get(hashing::hash_t hashed_kmer) const {
        return S.get_count(hashed_kmer);
    }

    /**
     * @Synopsis  Number of unique k-mers in the storage.
     *
     * @Returns   Number of unique k-mers.
     */
    uint64_t n_unique() const {
        return S.n_unique_kmers();
    }

    /**
     * @Synopsis  Number of occupied buckets in the storage.
     *
     * @Returns   Occupied buckets.
     */
    uint64_t n_occupied() const {
        return S.n_occupied();
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
     * @Synopsis  Given a root k-mers and the shift bases to its left,
     *            build the k-mer strings that could be its left neighbors.
     *
     * @Param nodes The shift_t objects containing the prefix bases and neighbor hashes.
     * @Param root The root k-mer.
     *
     * @Returns   kmer_t objects with the left k-mers.
     */
    std::vector<hashing::kmer_t> build_left_kmers(const std::vector<hashing::shift_t>& nodes,
                                                  const std::string& root) {
        std::vector<hashing::kmer_t> kmers;
        auto _prefix = prefix(root);
        for (auto neighbor : nodes) {
            kmers.push_back(hashing::kmer_t(neighbor.hash,
                                   neighbor.symbol
                                   + _prefix));
        }
        return kmers;
    }

    /**
     * @Synopsis  Given a root k-mers and the shift bases to its right,
     *            build the k-mer strings that could be its right neighbors.
     *
     * @Param nodes The shift_t obvjects containing the suffix bases and neighbor hashes.
     * @Param root The root k-mer.
     *
     * @Returns   kmer_t objects with the right k-mers.
     */
    std::vector<hashing::kmer_t> build_right_kmers(const std::vector<hashing::shift_t>& nodes,
                                                   const std::string& root) {
        std::vector<hashing::kmer_t> kmers;
        auto _suffix = suffix(root);
        for (auto neighbor : nodes) {
            kmers.push_back(hashing::kmer_t(neighbor.hash,
                                   _suffix
                                   + neighbor.symbol));
        }

        return kmers;
    }

    /**
     * @Synopsis  Finds the left-neighbors (in-neighbors) of root in the graph.
     *
     * @Param root The root k-mer.
     *
     * @Returns   shift_t objects with the prefix bases and hashes.
     */
    std::vector<hashing::shift_t> left_neighbors(const std::string& root) {
        auto assembler = this->get_assembler();
        assembler->set_cursor(root);
        return assembler->filter_nodes(assembler->gather_left());
    }

    /**
     * @Synopsis  Find the left-neighbors (in-neighbors) of root in the graph.
     *
     * @Param root The root k-mer.
     *
     * @Returns   kmer_t objects with the complete k-mers and hashes.
     */
    std::vector<hashing::kmer_t> left_neighbor_kmers(const std::string& root) {
        auto filtered = left_neighbors(root);
        return build_left_kmers(filtered, root);
    }

    /**
     * @Synopsis  Finds the right-neighbors (out-neighbors) of root in the graph.
     *
     * @Param root The root k-mer.
     *
     * @Returns   shift_t objects with the suffix bases and hashes.
     */
    std::vector<hashing::shift_t> right_neighbors(const std::string& root) {
        auto assembler = this->get_assembler();
        assembler->set_cursor(root);
        return assembler->filter_nodes(assembler->gather_right());
    }

    /**
     * @Synopsis  Finds the right-neighbors (out-neighbors) of root in the graph.
     *
     * @Param root The root k-mer.
     *
     * @Returns   kmer_t objects with the complete k-mers and hashes.
     */
    std::vector<hashing::kmer_t> right_neighbor_kmers(const std::string& root) {
        auto filtered = right_neighbors(root);
        return build_right_kmers(filtered, root);
    }

    /**
     * @Synopsis  Finds all neighbors of root in the graph.
     *
     * @Param root The root k-mer.
     *
     * @Returns   A pair of shift_t vectors; first contains the left and second the right neighbors.
     */
    std::pair<std::vector<hashing::shift_t>,
              std::vector<hashing::shift_t>> neighbors(const std::string& root) {
        auto assembler = this->get_assembler();
        assembler->set_cursor(root);
        auto lfiltered = assembler->filter_nodes(assembler->gather_left());
        auto rfiltered = assembler->filter_nodes(assembler->gather_right());
        return std::make_pair(lfiltered, rfiltered);
    }

    /**
     * @Synopsis  Finds all neighbors of root in the graph.
     *
     * @Param root The root k-mer.
     *
     * @Returns   A pair of kmer_t vectors; first contains the left and second the right neighbors.
     */
    std::pair<std::vector<hashing::kmer_t>,
              std::vector<hashing::kmer_t>> neighbor_kmers(const std::string& root) {
        auto filtered = neighbors(root);
        return std::make_pair(build_left_kmers(filtered.first, root),
                              build_right_kmers(filtered.second, root));
    }

    /**
     * @Synopsis  The address of the raw underling table, if available.
     *
     * @Returns   
     */
    uint8_t ** get_raw() const {
        return S.get_raw_tables();
    }

    /**
     * @Synopsis  If the StorageType is probabilistic, estimate its false positive rate.
     *
     * @Returns   The false-positive rate; 0 if an exact structure.
     */
    double estimated_fp() {
        double fp = n_occupied() / sizes[0];
        fp = pow(fp, N);
        return fp;
    }

    uint64_t add_sequence(const std::string& sequence,
                          std::vector<hashing::hash_t>& kmer_hashes,
                          std::vector<bool>& is_new) {
        hashing::KmerIterator<HashShifter> iter(sequence, _K);

        uint64_t n_consumed = 0;
        size_t pos = 0;
        bool kmer_consumed;
        while(!iter.done()) {
            hashing::hash_t h = iter.next();
            kmer_consumed = add(h);
            kmer_hashes.push_back(h);
            is_new.push_back(kmer_consumed);
            n_consumed += kmer_consumed;
            ++pos;
        }

        return n_consumed;
    }

    uint64_t add_sequence(const std::string& sequence,
                          std::set<hashing::hash_t>& new_kmers) {
        hashing::KmerIterator<HashShifter> iter(sequence, _K);

        uint64_t n_consumed = 0;
        size_t pos = 0;
        bool is_new;
        while(!iter.done()) {
            hashing::hash_t h = iter.next();
            is_new = add(h);
            if (is_new) {
                new_kmers.insert(h);
            }
            n_consumed += is_new;
            ++pos;
        }

        return n_consumed;
    }

    uint64_t add_sequence(const std::string& sequence) {
        hashing::KmerIterator<HashShifter> iter(sequence, _K);

        uint64_t n_consumed = 0;
        while(!iter.done()) {
            hashing::hash_t h = iter.next();
            n_consumed += add(h);
        }

        return n_consumed;
    }

    std::vector<storage::count_t> get_counts(const std::string& sequence) {
        hashing::KmerIterator<HashShifter> iter(sequence, _K);
        std::vector<storage::count_t> counts(sequence.length() - _K + 1);

        size_t pos = 0;
        while(!iter.done()) {
            hashing::hash_t h = iter.next();
            counts[pos] = get(h);
            ++pos;
        }

        return counts;
    }

    void get_counts(const std::string& sequence,
                    std::vector<storage::count_t>& counts,
                    std::vector<hashing::hash_t>& hashes,
                    std::set<hashing::hash_t>& new_hashes) {

        hashing::KmerIterator<HashShifter> iter(sequence, _K);

        while(!iter.done()) {
            hashing::hash_t h = iter.next();
            storage::count_t result = get(h);
            counts.push_back(result);
            hashes.push_back(h);
            if (result == 0) {
                new_hashes.insert(h);
            }
        }
    }

    void save(std::string filename) {
        S.save(filename, _K);
    }

    void load(std::string filename) {
        uint16_t ksize = _K;
        S.load(filename, ksize);
    }

    void reset() {
        S.reset();
    }

    std::shared_ptr<hashing::KmerIterator<HashShifter>> get_hash_iter(const std::string& sequence) {
        return std::make_shared<hashing::KmerIterator<HashShifter>>(sequence, _K);
    }

    std::shared_ptr<assembler_type> get_assembler() {
        auto ptr = this->shared_from_this();
        return std::make_shared<assembler_type>(ptr);
    }

};


}

#endif
