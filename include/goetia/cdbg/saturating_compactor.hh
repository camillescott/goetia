/**
 * (c) Camille Scott, 2019
 * File   : saturating_compactor.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.10.2019
 */

#ifndef GOETIA_SATURATING_CPTOR
#define GOETIA_SATURATING_CPTOR

#include "goetia/processors.hh"
#include "goetia/parsing/readers.hh"
#include "goetia/sequences/exceptions.hh"

namespace goetia {
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
                  uint64_t interval = metrics::IntervalCounter::DEFAULT_INTERVAL)
            : BaseType(interval),
              compactor(compactor),
              signature(signature)
        {

        }

        uint64_t process_sequence(const parsing::Record& read) {
            try {
                compactor->insert_sequence(read.sequence);
                signature->insert_sequence(read.sequence);
            } catch (InvalidCharacterException &e) {
                std::cerr << "WARNING: Bad sequence encountered at "
                          << this->_n_reads << ": "
                          << read.sequence << ", exception was "
                          << e.what() << std::endl;
                return 0;
            } catch (SequenceLengthException &e) {
                std::cerr << "NOTE: Skipped sequence that was too short: read "
                          << this->_n_reads << " with sequence "
                          << read.sequence 
                          << std::endl;
                return 0;
            } catch (std::exception &e) {
                std::cerr << "ERROR: Exception thrown at " << this->_n_reads 
                          << " with msg: " << e.what()
                          <<  std::endl;
                throw e;
            }

            return read.sequence.length() - compactor->K + 1;
        }

        void report() {

        }

        static std::shared_ptr<Processor> build(std::shared_ptr<compactor_type> compactor,
                                                std::shared_ptr<signature_type> signature,
                                                uint64_t interval = metrics::IntervalCounter::DEFAULT_INTERVAL) {
            return std::make_shared<Processor>(compactor,
                                               signature,
                                               interval);
        }

    };
};

}
}

#endif
