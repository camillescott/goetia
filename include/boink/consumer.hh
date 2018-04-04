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

template <class GraphType,
          class ParserType = FastxReader>
class FileConsumer {

    GraphType& graph;

public:

    FileConsumer(GraphType& graph) 
        : graph(graph) {
    }

    void consume(string const &filename,
                  uint64_t &n_reads,
                  uint64_t &n_consumed) {
        ReadParserPtr<ParserType> parser = get_parser<ParserType>(filename);
        consume(parser, n_reads, n_consumed);
    }

    void consume(read_parsers::ReadParserPtr<ParserType>& parser,
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
            uint64_t this_n_consumed = graph.add_sequence(read.cleaned_seq);

            __sync_add_and_fetch( &n_consumed, this_n_consumed );
            __sync_add_and_fetch( &n_reads, 1 );
        } 

    }

};

} //namespace boink
#endif
