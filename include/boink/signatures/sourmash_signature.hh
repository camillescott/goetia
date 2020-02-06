/**
 * (c) Camille Scott, 2019
 * File   : sourmash_signature.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 14.10.2019
 */

#ifndef BOINK_SOURMASH_SIG_HH
#define BOINK_SOURMASH_SIG_HH

#include <memory>

#include "boink/processors.hh"
#include "boink/parsing/readers.hh"
#include "boink/signatures/sourmash/kmer_min_hash.hh"


namespace boink {
namespace signatures {

struct SourmashSignature {

   class Signature: public sourmash::CKmerMinHash {
    
    public:

        const uint16_t K;

        Signature(unsigned int n, // n hashes
                  unsigned int K, // k-mer size
                  bool is_protein,
                  bool dayhoff,
                  bool hp,
                  uint32_t seed,
                  uint64_t max_hash)
            : sourmash::CKmerMinHash(n, K, is_protein, dayhoff, hp, seed, max_hash),
              K(K) 
        {
        }

        Signature(CKmerMinHash* minhash)
            : sourmash::CKmerMinHash(minhash->_get_ptr()),
              K(K)
        {
        }

        static std::shared_ptr<Signature> build(unsigned int n, // n hashes
                                                unsigned int K, // k-mer size
                                                bool is_protein,
                                                bool dayhoff,
                                                bool hp,
                                                uint32_t seed,
                                                uint64_t max_hash) {
            return std::make_shared<Signature>(n, K, is_protein,
                                               dayhoff, hp,
                                               seed, max_hash);
        }

        size_t insert_sequence(const std::string& sequence) {
            this->add_sequence(sequence.c_str());
            return 1;
        }
    };


    using Processor = InserterProcessor<Signature, parsing::FastxParser<DNAN_SIMPLE>>;


};


}
}

#endif
