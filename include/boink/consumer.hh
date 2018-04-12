/* consumer.hh -- dBG consumers
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef CONSUMER_HH
#define CONSUMER_HH

#include <memory>
#include <iostream>

#include "oxli/read_parsers.hh"

using namespace oxli;
using namespace oxli:: read_parsers;

#define DEFAULT_OUTPUT_INTERVAL 10000

namespace boink {

template <class Derived,
          class ParserType = FastxReader>
class FileProcessor {

public:

    void process(string const &filename,
                  uint64_t &n_reads,
                  uint64_t &n_consumed) {
        ReadParserPtr<ParserType> parser = get_parser<ParserType>(filename);
        process(parser, n_reads, n_consumed);
    }

    void process(string const &filename,
                  uint64_t &n_reads,
                  uint64_t &n_consumed,
                  uint32_t output_interval) {
        ReadParserPtr<ParserType> parser = get_parser<ParserType>(filename);
        process(parser, n_reads, n_consumed, output_interval);
    }

    void process(read_parsers::ReadParserPtr<ParserType>& parser,
                 uint64_t &n_reads,
                 uint64_t &n_consumed) {
        Read read;

        // Iterate through the reads and consume their k-mers.
        while (!parser->is_complete( )) {
            try {
                read = parser->get_next_read( );
            } catch (NoMoreReadsAvailable) {
                break;
            }

            read.set_clean_seq();
            uint64_t this_n_consumed = derived().process_sequence(read.cleaned_seq,
                                                                  n_reads);
                
                //graph.add_sequence(read.cleaned_seq);

            __sync_add_and_fetch( &n_consumed, this_n_consumed );
            __sync_add_and_fetch( &n_reads, 1 );
        }
    }

    void process(read_parsers::ReadParserPtr<ParserType>& parser,
                 uint64_t &n_reads,
                 uint64_t &n_consumed,
                 uint32_t output_interval) {
        Read read;

        // Iterate through the reads and consume their k-mers.
        uint32_t output_counter = 0;
        while (!parser->is_complete( )) {
            try {
                read = parser->get_next_read( );
            } catch (NoMoreReadsAvailable) {
                break;
            }

            read.set_clean_seq();
            uint64_t this_n_consumed = derived().process_sequence(read,
                                                                  n_reads);
                
            __sync_add_and_fetch( &n_consumed, this_n_consumed );
            __sync_add_and_fetch( &n_reads, 1 );
            __sync_add_and_fetch( &output_counter, 1 );

            if (output_counter == output_interval) {
                output_counter = 0;
                std::cerr << "processed " << n_reads << " reads, "
                          << n_consumed << " unique k-mers." << std::endl;
            }    
        }
    }

private:

    friend Derived;

    Derived& derived() {
        return *static_cast<Derived*>(this);
    }

    const Derived& derived() const {
        return *static_cast<const Derived*>(this);
    }

};


template <class GraphType,
          class ParserType = FastxReader>
class FileConsumer : public FileProcessor<FileConsumer<GraphType, ParserType>,
                                          ParserType> {

protected:

    GraphType * graph;

public:

    FileConsumer(GraphType * graph)
        : graph(graph) {

    }

    uint64_t process_sequence(const Read& read,
                              const uint64_t read_n) {
        return graph->add_sequence(read.cleaned_seq);
    }

};


template <class GraphType,
          class ParserType = FastxReader>
class DecisionNodeProcessor : public FileProcessor<DecisionNodeProcessor<GraphType, ParserType>,
                                                   ParserType> {

protected:

    StreamingCompactor<GraphType> * compactor;
    GraphType * graph;
    std::string _output_filename;
    std::ofstream _output_stream;

public:

    DecisionNodeProcessor(StreamingCompactor<GraphType> * compactor,
                          std::string& output_filename)
        : compactor(compactor),
          graph(compactor->dbg),
          _output_filename(output_filename),
          _output_stream(_output_filename.c_str()) {

        _output_stream << "read_n, l_degree, r_degree, position, hash" << std::endl;

    }

    ~DecisionNodeProcessor() {
        _output_stream.close();
    }

    uint64_t process_sequence(const Read& read,
                              const uint64_t read_n) {
        uint64_t n_new = graph->add_sequence(read.cleaned_seq);
        if (n_new > 0) {
            std::vector<uint32_t> decision_positions;
            HashVector decision_hashes;
            std::vector<NeighborBundle> decision_neighbors;
            std::set<hash_t> new_kmers;

            compactor->find_decision_nodes(read.cleaned_seq,
                                           decision_positions,
                                           decision_hashes,
                                           decision_neighbors);

            for (size_t i=0; i<decision_positions.size(); ++i) {
                _output_stream << read_n << ", "
                               << decision_neighbors[i].first.size() << ", "
                               << decision_neighbors[i].second.size() << ", "
                               << decision_positions[i] << ", "
                               << decision_hashes[i] << std::endl;
            }
        }

        return n_new;
    }



};


template <class ShifterType,
          class ParserType = FastxReader>
class MinimizerProcessor : public FileProcessor<MinimizerProcessor<ShifterType, ParserType>,
                                                ParserType> {

protected:

    WKMinimizer<ShifterType> M;
    std::string _output_filename;
    std::ofstream _output_stream;

public:

    MinimizerProcessor(int32_t window_size,
                       uint16_t K,
                       const std::string& output_filename)
        : M(window_size, K),
          _output_filename(output_filename),
          _output_stream(_output_filename.c_str()) {
    
    }

    ~MinimizerProcessor() {
        _output_stream.close();
    }

    uint64_t process_sequence(const Read& read,
                              const uint64_t read_n) {
        std::vector<typename WKMinimizer<ShifterType>::value_type> minimizers;
        minimizers = M.get_minimizers(read.cleaned_seq);

        for (auto min : minimizers) {
            _output_stream << read_n << ","
                           << min.second << ","
                           << min.first << ","
                           << read.cleaned_seq.substr(min.second, M.K())
                           << std::endl;
        }
    }
};

} //namespace boink
#endif
