/* diginorm.hh -- 
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_DIGINORM_HH
#define BOINK_DIGINORM_HH

#include "boink/processors.hh"
#include "boink/dbg.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/hashing/exceptions.hh"

namespace boink {
namespace normalization {


template <class ShifterType>
bool median_count_at_least(const std::string&          sequence,
                           unsigned int                cutoff,
                           dBG<storage::ByteStorage,
                               ShifterType>          * counts) {

    auto kmers = counts->get_hash_iter(sequence);
    unsigned int min_req = 0.5 + float(sequence.size() - counts->K() + 1) / 2;
    unsigned int num_cutoff_kmers = 0;

    // first loop:
    // accumulate at least min_req worth of counts before checking to see
    // if we have enough high-abundance k-mers to indicate success.
    for (unsigned int i = 0; i < min_req; ++i) {
        hash_t kmer = kmers->next();
        if (counts->query(kmer) >= cutoff) {
            ++num_cutoff_kmers;
        }
    }


    // second loop: now check to see if we pass the threshold for each k-mer.
    if (num_cutoff_kmers >= min_req) {
        return true;
    }
    while(!kmers->done()) {
        hash_t kmer = kmers->next();
        if (counts->query(kmer) >= cutoff) {
            ++num_cutoff_kmers;
            if (num_cutoff_kmers >= min_req) {
                return true;
            }
        }
    }
    return false;
}


template <class GraphType,
          class ParserType = parsing::FastxReader>
class NormalizingCompactor: 
    public FileProcessor<NormalizingCompactor<GraphType, ParserType>,
                         ParserType> { //template class names like modern art

protected:

    std::shared_ptr<cdbg::StreamingCompactor<GraphType>>   compactor;
    std::shared_ptr<GraphType>                             graph;

    std::unique_ptr<dBG<storage::ByteStorage,
                        typename GraphType::shifter_type>> counts;
    unsigned int                                           cutoff;
    size_t                                                 n_seq_updates;

    typedef FileProcessor<NormalizingCompactor<GraphType, ParserType>,
                          ParserType> Base;

public:

    using Base::process_sequence;
    using events::EventNotifier::register_listener;
    
    NormalizingCompactor(shared_ptr<cdbg::StreamingCompactor<GraphType>> compactor,
                         unsigned int cutoff,
                         uint64_t fine_interval=DEFAULT_FINE_INTERVAL,
                         uint64_t medium_interval=DEFAULT_MEDIUM_INTERVAL,
                         uint64_t coarse_interval=DEFAULT_COARSE_INTERVAL)
        : Base(fine_interval, medium_interval, coarse_interval),
          compactor(compactor),
          graph(compactor->dbg),
          cutoff(cutoff),
          n_seq_updates(0)
    {
        counts = std::make_unique<dBG<storage::ByteStorage,
                                      typename GraphType::shifter_type>>(graph->K(),
                                                                         100000000,
                                                                         4);
    }

    void process_sequence(const parsing::Read& read) {

        if (median_count_at_least(read.cleaned_seq, cutoff, counts.get())) {
            return;
        }

        counts->insert_sequence(read.cleaned_seq);

        try {
            compactor->update_sequence(read.cleaned_seq);
        } catch (hashing::InvalidCharacterException &e) {
            std::cerr << "WARNING: Bad sequence encountered at "
                      << this->_n_reads << ": "
                      << read.cleaned_seq << ", exception was "
                      << e.what() << std::endl;
            return;
        } catch (hashing::SequenceLengthException &e) {
            std::cerr << "NOTE: Skipped sequence that was too short: read "
                      << this->_n_reads << " with sequence "
                      << read.cleaned_seq 
                      << std::endl;
            return;
        } catch (std::exception &e) {
            std::cerr << "ERROR: Exception thrown at " << this->_n_reads 
                      << " with msg: " << e.what()
                      <<  std::endl;
            throw e;
        }

        ++n_seq_updates;
    }

    void report() {
        std::cerr << "\t" << n_seq_updates << " used for cDBG updates." << std::endl;
    }

};

}
}



#endif
