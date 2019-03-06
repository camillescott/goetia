/* ukhs_signature.hh
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_UKHS_SIGNATURE_HH
#define BOINK_UKHS_SIGNATURE_HH

#include <algorithm>
#include <memory>
#include <vector>

#include "boink/boink.hh"
#include "boink/hashing/alphabets.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/hashshifter.hh"
#include "boink/kmers/kmerclient.hh"
#include "boink/hashing/ukhs.hh"

namespace boink {
namespace signatures {

using boink::hashing::hash_t;
using boink::hashing::PartitionedHash;


class IncompatibleSignature : public BoinkException {
public:
    explicit IncompatibleSignature(const std::string& msg = "Incompatible signatures.")
        : BoinkException(msg) { }
};


struct MinHash {
    // this based off and partially riffed from sourmash
    std::vector<hash_t> mins;
    const size_t max_size;
    uint64_t     n_inserts;

    explicit MinHash(size_t max_size) 
        : max_size(max_size),
          mins(max_size, ULLONG_MAX),
          n_inserts(0)
    {
    }

    bool insert(const hash_t hash) {
        ++n_inserts;
        if (mins.size() == 0) {
            mins.push_back(hash);
            return true;
        } else if (mins.back() > hash or mins.size() < max_size) {
            auto pos = std::lower_bound(std::begin(mins), std::end(mins), hash);

            if (pos == mins.cend()) {
                mins.push_back(hash);
                return true;
            } else if (*pos != hash) {
                mins.insert(pos, hash);
            
                if (max_size and mins.size() > max_size) {
                    mins.pop_back();
                }

                return true;
            }
        }

        return false;
        
    }
};


class UKHSSignature : public kmers::KmerClient {

protected:

    std::shared_ptr<hashing::UKHS>        ukhs;
    hashing::DefaultUKHSShifter           partitioner;

    std::vector<std::unique_ptr<MinHash>> sig_buckets;

public:

    const uint64_t                        n_buckets;
    const uint16_t                        bucket_K;
    const uint64_t                        bucket_size;

    uint64_t                              n_accepted;
    uint64_t                              n_rejected;
    uint64_t                              n_kmers;

    explicit UKHSSignature(uint16_t K,
                           uint16_t bucket_K,
                           uint64_t bucket_size,
                           std::shared_ptr<hashing::UKHS> ukhs)
        : KmerClient  (K),
          ukhs        (ukhs),
          partitioner (K, bucket_K, ukhs),
          n_buckets   (partitioner.n_ukhs_hashes()),
          bucket_K    (bucket_K),
          bucket_size (bucket_size),
          n_accepted  (0),
          n_rejected  (0),
          n_kmers     (0)
    {
        for (size_t i = 0; i < n_buckets; ++i) {
            sig_buckets.push_back(
                std::make_unique<MinHash>(bucket_size)
            );
            for (size_t j = 0; j < bucket_size; ++j) {
                sig_buckets.back()->insert(ULLONG_MAX);
            }
        }
    }

    inline void insert(const std::string& kmer) {
        partitioner.set_cursor(kmer);
        partitioner.reset_unikmers();
        ++n_kmers;

        insert(partitioner.get(), partitioner.get_front_partition());
    }

    inline void insert(const PartitionedHash& ph) {
        insert(ph.first, ph.second);
    }

    inline void insert(const hash_t hash, const uint64_t bucket_id) {
        if (this->query_bucket(bucket_id)->insert(hash)) {
            ++n_accepted;
        } else {
            ++n_rejected;
        }
    }

    inline void insert_sequence(const std::string& sequence) {
        hashing::KmerIterator<hashing::DefaultUKHSShifter> iter(sequence, &partitioner);

        PartitionedHash h          = iter.next();
        uint64_t        cur_bid    = h.second;
        MinHash *       cur_bucket = this->query_bucket(cur_bid);
        if (cur_bucket->insert(h.first)) {
            ++n_accepted;
        } else {
            ++n_rejected;
        }
        ++n_kmers;

        while(!iter.done()) {
            h = iter.next();
            if (h.second != cur_bid) {
                cur_bid = h.second;
                cur_bucket = this->query_bucket(cur_bid);
            }
            if (cur_bucket->insert(h.first)) {
                n_accepted++;
            } else {
                n_rejected++;
            }
            ++n_kmers;
        }
    }

    std::vector<std::vector<hash_t>> get_signature() {
        std::vector<std::vector<hash_t>> signature;
        for (auto& bucket : sig_buckets) {
            signature.push_back(std::vector<hash_t>(bucket->mins));
        }

        return signature;
    }

    std::vector<uint64_t> get_bucket_n_inserts() {
        std::vector<uint64_t> result;
        for (auto& bucket : sig_buckets) {
            result.push_back(bucket->n_inserts);
        }
        return result;
    }

    uint64_t get_n_accepted() const {
        return n_accepted;
    }

    uint64_t get_n_rejected() const {
        return n_rejected;
    }

    uint64_t get_n_kmers() const {
        return n_kmers;
    }

    std::set<hash_t> intersection(UKHSSignature * other) {
        if (other->n_buckets != this->n_buckets or
            other->bucket_K  != this->bucket_K or
            other->K()       != this->K() or
            other->n_buckets != this->n_buckets) {
            
            throw IncompatibleSignature("Error: Signatures not compatible");
        }
    }

protected:

    MinHash * query_bucket(uint64_t bucket_id) {
        if (bucket_id < n_buckets) {
            return sig_buckets[bucket_id].get();
        } else {
            throw BoinkException("Invalid UKHS bucket: " + std::to_string(bucket_id));
        }
    }
          
};

}
}
#endif
