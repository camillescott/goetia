/**
 * (c) Camille Scott, 2019
 * File   : ukhs_signature.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
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
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/hashextender.hh"
#include "boink/event_types.hh"
#include "boink/reporting/reporters.hh"
#include "boink/processors.hh"
#include "boink/kmers/kmerclient.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/hashing/ukhshashshifter.hh"
#include "boink/hashing/canonical.hh"

#include "boink/pdbg.hh"

namespace boink {
namespace signatures {


class IncompatibleSignature : public BoinkException {
public:
    explicit IncompatibleSignature(const std::string& msg = "Incompatible signatures.")
        : BoinkException(msg) { }
};


template <class StorageType, class BaseShifterType>
struct UnikmerSignature {

    typedef StorageType                              storage_type;
    typedef hashing::UnikmerShifter<BaseShifterType> shifter_type;
    typedef typename shifter_type::ukhs_type         ukhs_type;
    typedef PdBG<StorageType, BaseShifterType>       pdbg_type;

    _boink_model_typedefs_from_shiftertype(shifter_type)

    class Signature : public kmers::KmerClient {

    protected:

        std::shared_ptr<ukhs_type> ukhs_map;

        std::shared_ptr<pdbg_type> signature;

    public:

        const uint16_t                        bucket_K;

        template <typename... Args>
        explicit Signature(uint16_t K,
                           uint16_t bucket_K,
                           std::shared_ptr<ukhs_type> ukhs_map,
                           Args&&... args)
            : KmerClient  (K),
              ukhs_map    (ukhs_map),
              bucket_K    (bucket_K)
        {
            signature = std::make_shared<pdbg_type>(K,
                                                    bucket_K,
                                                    ukhs_map,
                                                    std::forward<Args>(args)...);
        }

        template<typename...Args>
        static std::shared_ptr<Signature> build(uint16_t K,
                                                uint16_t bucket_K,
                                                std::shared_ptr<ukhs_type> ukhs_map,
                                                Args&&... args) {
            return std::make_shared<Signature>(K, bucket_K, ukhs_map, std::forward<Args>(args)...);
        }

        template<typename U = StorageType>
        static std::shared_ptr<Signature> build(uint16_t K,
                                                uint16_t bucket_K,
                                                std::shared_ptr<ukhs_type> ukhs_map,
                                                typename std::enable_if_t<std::is_same<U, boink::storage::SparseppSetStorage>::value, U*> = 0) {
            return std::make_shared<Signature>(K, bucket_K, ukhs_map);
        }

        inline void insert(const std::string& kmer) {
            signature->insert(kmer);
        }

        inline size_t insert_sequence(const std::string& sequence) {
            return signature->insert_sequence(sequence);
        }

        size_t get_size() const {
            return signature->n_partitions();
        }

        std::vector<size_t> get_signature() {
            return signature->get_partition_counts();
        }

        uint64_t get_n_kmers() const {
            size_t n_kmers = 0;
            for (const auto& c : signature->get_partition_counts()) {
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

    typedef UnikmerSignature<StorageType, BaseShifterType>::Signature signature_type;

    using Processor = InserterProcessor<Signature>; 


    class Reporter : public reporting::SingleFileReporter {

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

};

}
}


#endif
