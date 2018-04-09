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

    uint64_t process_sequence(const std::string& sequence,
                              const uint64_t read_n) {
        return graph->add_sequence(sequence);
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

        _output_stream << "read_n, l_degree, r_degree, hash, position" << std::endl;

    }

    ~DecisionNodeProcessor() {
        _output_stream.close();
    }

    uint64_t process_sequence(const std::string& sequence,
                              const uint64_t read_n) {
        uint64_t n_new = graph->add_sequence(sequence);
        if (n_new > 0) {
            std::vector<uint32_t> decision_positions;
            HashVector decision_hashes;
            std::vector<NeighborBundle> decision_neighbors;

            compactor->find_decision_nodes(sequence,
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

} //namespace boink
#endif
