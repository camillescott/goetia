/**
 * (c) Camille Scott, 2019
 * File   : ukhs_signature.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */
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
#include <fstream>
#include <iostream>
#include <string>
#include <memory>
#include <set>
#include <vector>

#include "boink/boink.hh"
#include "boink/hashing/alphabets.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/event_types.hh"
#include "boink/reporting/reporters.hh"
#include "boink/processors.hh"
#include "boink/kmers/kmerclient.hh"
#include "boink/hashing/ukhs.hh"

#include "boink/pdbg.hh"

namespace boink {
namespace signatures {


class IncompatibleSignature : public BoinkException {
public:
    explicit IncompatibleSignature(const std::string& msg = "Incompatible signatures.")
        : BoinkException(msg) { }
};


template <class StorageType>
struct UnikmerSignature {

    typedef StorageType                      storage_type;
    typedef hashing::UKHS::LazyShifter       shifter_type;
    typedef typename shifter_type::hash_type hash_type;
    typedef typename shifter_type::kmer_type kmer_type;
    typedef typename shifter_type::shift_type shift_type;

    class Signature : public kmers::KmerClient {

    protected:

        std::shared_ptr<hashing::UKHS::Map> ukhs_map;

        std::shared_ptr<PdBG<StorageType>> signature;

    public:

        const uint16_t                        bucket_K;

        template <typename... Args>
        explicit Signature(uint16_t K,
                           uint16_t bucket_K,
                           std::shared_ptr<hashing::UKHS::Map> ukhs_map,
                           Args&&... args)
            : KmerClient  (K),
              ukhs_map    (ukhs_map),
              bucket_K    (bucket_K)
        {
            signature = std::make_shared<PdBG<StorageType>>(K,
                                                            bucket_K,
                                                            ukhs_map,
                                                            std::forward<Args>(args)...);
        }

        template<typename...Args>
        static std::shared_ptr<Signature> build(uint16_t K,
                                                uint16_t bucket_K,
                                                std::shared_ptr<hashing::UKHS::Map> ukhs_map,
                                                Args&&... args) {
            return std::make_shared<Signature>(K, bucket_K, ukhs_map, std::forward<Args>(args)...);
        }

        template<typename U = StorageType>
        static std::shared_ptr<Signature> build(uint16_t K,
                                                uint16_t bucket_K,
                                                std::shared_ptr<hashing::UKHS::Map> ukhs_map,
                                                typename std::enable_if_t<std::is_same<U, boink::storage::SparseppSetStorage<>>::value, U*> = 0) {
            return std::make_shared<Signature>(K, bucket_K, ukhs_map);
        }

        inline void insert(const std::string& kmer) {
            signature->insert(kmer);
        }

        inline void insert_sequence(const std::string& sequence) {
            signature->insert_sequence(sequence);
        }

        size_t get_size() const {
            return signature->n_partitions();
        }

        std::vector<size_t> get_signature() {
            return signature->get_partition_counts();
        }

        uint64_t get_n_kmers() const {
            size_t n_kmers = 0;
            for (auto& c : signature->get_partition_counts()) {
                n_kmers += c;
            }
            return n_kmers;
        }

        std::set<hash_type> intersection(Signature * other) {
            if (other->get_size() != this->get_size() or
                other->bucket_K   != this->bucket_K or
                other->K()        != this->K()) {
                
                throw IncompatibleSignature("Error: Signatures not compatible");
            }
        }
    };

    /*
    class Reporter : public reporting::SingleFileReporter {
    public:

        Reporter(const std::string& filename)
            : SingleFileReporter(filename, "Test")
        {
        }

        virtual void handle_msg(std::shared_ptr<events::Event> event) {

        }

    };
    */    


    class Processor : public FileProcessor<Processor,
                                           parsing::FastxReader> {
    protected:

        std::shared_ptr<Signature> signature;

        typedef FileProcessor<Processor, parsing::FastxReader> Base;

    public:

        using Base::process_sequence;

        Processor(std::shared_ptr<Signature> signature,
                  uint64_t fine_interval   ,
                  uint64_t medium_interval ,
                  uint64_t coarse_interval )
            : Base(fine_interval, medium_interval, coarse_interval),
              signature(signature)
        {
        }

        static std::shared_ptr<Processor> build(std::shared_ptr<Signature> signature,
                                                uint64_t fine_interval   = DEFAULT_INTERVALS::FINE,
                                                uint64_t medium_interval = DEFAULT_INTERVALS::MEDIUM,
                                                uint64_t coarse_interval = DEFAULT_INTERVALS::COARSE) {
            return std::make_shared<Processor>(signature, fine_interval, medium_interval, coarse_interval);
        }

        void process_sequence(const parsing::Read& read) {
            signature->insert_sequence(read.cleaned_seq);
        }

        void report() {
        }
    };

/*
    class Reporter : public reporting::SingleFileReporter; {

    private:

        std::shared_ptr<Signature> signature;

    public:

        Reporter(std::shared_ptr<Signature> signature,
                 const std::string&         filename);

        
        static std::shared_ptr<Reporter> build(std::shared_ptr<Signature> signature,
                                               const std::string&         filename) {
            return std::make_shared<Reporter>(signature, filename);
        }
        
        
        virtual void handle_msg(std::shared_ptr<events::Event> event);
    
    };
*/

};

}
}


#endif
