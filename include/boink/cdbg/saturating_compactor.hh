/**
 * (c) Camille Scott, 2019
 * File   : saturating_compactor.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.10.2019
 */

#ifndef BOINK_SATURATING_CPTOR
#define BOINK_SATURATING_CPTOR

#include "boink/processors.hh"
#include "boink/parsing/readers.hh"
#include "boink/sequences/exceptions.hh"

namespace boink {
namespace cdbg {

template <class CompactorType,
          class SignatureType>
struct SaturatingCompactor {

    typedef CompactorType                       compactor_type;
    typedef SignatureType                       signature_type;

    typedef typename compactor_type::hash_type  hash_type;
    typedef typename compactor_type::kmer_type  kmer_type;


    class Processor : public FileProcessor<Processor,
                                           parsing::FastxParser<>> {

    protected:

        typedef FileProcessor<Processor, parsing::FastxParser<>> BaseType;

    public:
        
        std::shared_ptr<compactor_type> compactor;
        std::shared_ptr<signature_type> signature;

        using BaseType::process_sequence;

        Processor(std::shared_ptr<compactor_type> compactor,
                  std::shared_ptr<signature_type> signature,
                  uint64_t fine_interval   = DEFAULT_INTERVALS::FINE,
                  uint64_t medium_interval = DEFAULT_INTERVALS::MEDIUM,
                  uint64_t coarse_interval = DEFAULT_INTERVALS::COARSE)
            : BaseType(fine_interval, medium_interval, coarse_interval),
              compactor(compactor),
              signature(signature)
        {

        }

        void process_sequence(const parsing::Record& read) {
            try {
                compactor->insert_sequence(read.sequence);
                signature->insert_sequence(read.sequence);
            } catch (InvalidCharacterException &e) {
                std::cerr << "WARNING: Bad sequence encountered at "
                          << this->_n_reads << ": "
                          << read.sequence << ", exception was "
                          << e.what() << std::endl;
                                              return;
            } catch (SequenceLengthException &e) {
                std::cerr << "NOTE: Skipped sequence that was too short: read "
                          << this->_n_reads << " with sequence "
                          << read.sequence 
                          << std::endl;
                return;
            } catch (std::exception &e) {
                std::cerr << "ERROR: Exception thrown at " << this->_n_reads 
                          << " with msg: " << e.what()
                          <<  std::endl;
                throw e;
            }
        }

        void report() {

        }

        static std::shared_ptr<Processor> build(std::shared_ptr<compactor_type> compactor,
                                                std::shared_ptr<signature_type> signature,
                                                uint64_t fine_interval   = DEFAULT_INTERVALS::FINE,
                                                uint64_t medium_interval = DEFAULT_INTERVALS::MEDIUM,
                                                uint64_t coarse_interval = DEFAULT_INTERVALS::COARSE) {
            return std::make_shared<Processor>(compactor,
                                               signature,
                                               fine_interval,
                                               medium_interval,
                                               coarse_interval);
        }

    };
};

}
}

#endif
