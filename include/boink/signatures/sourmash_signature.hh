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

#include "boink/kmers/kmerclient.hh"
#include "boink/processors.hh"
#include "boink/signatures/sourmash/kmer_min_hash.hh"


namespace boink {
namespace signatures {

struct SourmashSignature {

    class Signature: public KmerMinHash, 
                     public kmers::KmerClient {
    
    public:

        Signature(unsigned int n, // n hashes
                  unsigned int K, // k-mer size
                  bool is_protein,
                  uint32_t seed,
                  uint64_t max_hash)
            : KmerMinHash(n, K, is_protein, seed, max_hash),
              KmerClient(K) {
            
        }

        static std::shared_ptr<Signature> build(unsigned int n, // n hashes
                                                unsigned int K, // k-mer size
                                                bool is_protein,
                                                uint32_t seed,
                                                uint64_t max_hash) {
            return std::make_shared<Signature>(n, K, is_protein,
                                               seed, max_hash);
        }

        size_t insert_sequence(const std::string& sequence) {
            this->add_sequence(sequence.c_str());
            return 1;
        }
    };


    using Processor = InserterProcessor<Signature>;


};


}
}

#endif
