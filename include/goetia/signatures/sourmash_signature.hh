/**
 * (c) Camille Scott, 2019
 * File   : sourmash_signature.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 14.10.2019
 */

#ifndef GOETIA_SOURMASH_SIG_HH
#define GOETIA_SOURMASH_SIG_HH

#include <memory>

#include "goetia/parsing/readers.hh"
#include "goetia/processors.hh"
#include "goetia/sequences/alphabets.hh"
#include "goetia/signatures/sourmash/sourmash.hpp"


namespace goetia::signatures {

    struct SourmashSignature {

        class Signature : public sourmash::MinHash {

          public:
            const uint16_t K;

            Signature(unsigned int n,  // n hashes
                      unsigned int K,  // k-mer size
                      bool         is_protein,
                      bool         dayhoff,
                      bool         hp,
                      uint32_t     seed,
                      uint64_t     scaled)
                : sourmash::MinHash(n, K, is_protein, dayhoff, hp, seed, Signature::max_hash_from_scaled(scaled)), K(K) {
            }

            Signature(MinHash* minhash) : sourmash::MinHash(minhash->_get_ptr()), K(minhash->ksize()) {
            }

            static std::shared_ptr<Signature> build(unsigned int n,  // n hashes
                                                    unsigned int K,  // k-mer size
                                                    bool         is_protein,
                                                    bool         dayhoff,
                                                    bool         hp,
                                                    uint32_t     seed,
                                                    uint64_t     max_hash) {
                return std::make_shared<Signature>(n, K, is_protein, dayhoff, hp, seed, max_hash);
            }

            static uint64_t max_hash_from_scaled(uint64_t scaled) {
                if (scaled == 0) {
                    return 0;
                } else if (scaled == 1) {
                    return std::numeric_limits<uint64_t>::max();
                } else {
                    return static_cast<uint64_t>(static_cast<double>(std::numeric_limits<uint64_t>::max()) /
                                                 static_cast<double>(scaled));
                }
            }

            static uint64_t scaled_from_max_hash(uint64_t max_hash) {
                if (max_hash == 0) {
                    return 0;
                } else {
                    return std::numeric_limits<uint64_t>::max() / max_hash;
                }
            }

            size_t insert_sequence(const std::string& sequence) {
                this->add_sequence(sequence.c_str());
                return sequence.length() - K + 1;
            }
        };

        using Processor = InserterProcessor<Signature, parsing::FastxParser<DNAN_SIMPLE>>;
    };
}  // namespace goetia::signatures

#endif
