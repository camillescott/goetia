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

#include "goetia/goetia.hh"
#include "goetia/event_types.hh"
#include "goetia/reporting/reporters.hh"
#include "goetia/processors.hh"

#include "goetia/hashing/kmeriterator.hh"
#include "goetia/hashing/hashextender.hh"
#include "goetia/hashing/ukhs.hh"
#include "goetia/hashing/unikmershifter.hh"
#include "goetia/hashing/canonical.hh"
#include "goetia/storage/storage_types.hh"

#include "goetia/pdbg.hh"

namespace goetia {
namespace signatures {


class IncompatibleSignature : public BoinkException {
public:
    explicit IncompatibleSignature(const std::string& msg = "Incompatible signatures.")
        : BoinkException(msg) { }
};


template <class StorageType, class HashType>
struct UnikmerSignature {

    typedef StorageType                                 storage_type;
    typedef hashing::UnikmerLemirePolicy<HashType>      shift_policy;
    typedef hashing::HashShifter<shift_policy>          shifter_type;
    typedef typename shifter_type::ukhs_type            ukhs_type;
    typedef PdBG<StorageType, shifter_type>             pdbg_type;

    _goetia_model_typedefs_from_shiftertype(shifter_type)

    class Signature {

    protected:

        std::shared_ptr<ukhs_type> ukhs_map;
        std::shared_ptr<pdbg_type> signature;

    public:

        const uint16_t W;
        const uint16_t K;

        template <typename... Args>
        explicit Signature(uint16_t W,
                           uint16_t K,
                           std::shared_ptr<ukhs_type> ukhs_map,
                           Args&&... args)
            : W        (W),
              K        (K),
              ukhs_map (ukhs_map)
        {
            signature = std::make_shared<pdbg_type>(W,
                                                    K,
                                                    ukhs_map,
                                                    std::forward<Args>(args)...);
        }

        template<typename...Args>
        static std::shared_ptr<Signature> build(uint16_t W,
                                                uint16_t K,
                                                std::shared_ptr<ukhs_type> ukhs_map,
                                                Args&&... args) {
            return std::make_shared<Signature>(W, K, ukhs_map, std::forward<Args>(args)...);
        }

        template<typename U = StorageType>
        static std::shared_ptr<Signature> build(uint16_t W,
                                                uint16_t K,
                                                std::shared_ptr<ukhs_type> ukhs_map,
                                                typename std::enable_if_t<std::is_same<U, goetia::storage::SparseppSetStorage>::value, U*> = 0) {
            return std::make_shared<Signature>(W, K, ukhs_map);
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

        /*
        std::set<hash_type> intersection(Signature * other) {
            if (other->get_size() != this->get_size() or
                other->bucket_K   != this->bucket_K or
                other->K()        != this->K()) {
                
                throw IncompatibleSignature("Error: Signatures not compatible");
            }
        }
        */
    };

    typedef typename UnikmerSignature<StorageType, HashType>::Signature signature_type;

    using Processor = InserterProcessor<Signature>; 


    class Reporter : public reporting::SingleFileReporter {

    private:

        std::shared_ptr<Signature> signature;

    public:

        Reporter(std::shared_ptr<Signature> signature,
                 const std::string&         filename)
            : SingleFileReporter(filename, "UnikmerSignature::Reporter"),
              signature(signature)
        {
            _cerr(this->THREAD_NAME << " reporting at MEDIUM interval.");
            this->msg_type_whitelist.insert(events::MSG_TIME_INTERVAL);
        }

        
        static std::shared_ptr<Reporter> build(std::shared_ptr<Signature> signature,
                                               const std::string&         filename) {
            return std::make_shared<Reporter>(signature, filename);
        }
        
        
        virtual void handle_msg(std::shared_ptr<events::Event> event) {
             if (event->msg_type == events::MSG_TIME_INTERVAL) {
                auto _event = static_cast<events::TimeIntervalEvent*>(event.get());
                if (_event->level == events::TimeIntervalEvent::MEDIUM ||
                    _event->level == events::TimeIntervalEvent::END) {
                    
                    _output_stream << _event->t;
                    auto counts = signature->get_signature();
                    for (auto& count : counts) {
                        _output_stream << ", " << count;
                    }
                    _output_stream << std::endl;
                }
            }       
        }
    
    };

};

extern template class signatures::UnikmerSignature<storage::BitStorage, hashing::Hash<uint64_t>>;
extern template class signatures::UnikmerSignature<storage::BitStorage, hashing::Canonical<uint64_t>>;

extern template class signatures::UnikmerSignature<storage::SparseppSetStorage, hashing::Hash<uint64_t>>;
extern template class signatures::UnikmerSignature<storage::SparseppSetStorage, hashing::Canonical<uint64_t>>;

extern template class signatures::UnikmerSignature<storage::ByteStorage, hashing::Hash<uint64_t>>;
extern template class signatures::UnikmerSignature<storage::ByteStorage, hashing::Canonical<uint64_t>>;

extern template class signatures::UnikmerSignature<storage::NibbleStorage, hashing::Hash<uint64_t>>;
extern template class signatures::UnikmerSignature<storage::NibbleStorage, hashing::Canonical<uint64_t>>;

extern template class signatures::UnikmerSignature<storage::QFStorage, hashing::Hash<uint64_t>>;
extern template class signatures::UnikmerSignature<storage::QFStorage, hashing::Canonical<uint64_t>>;



}
}


#endif
