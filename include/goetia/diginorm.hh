/**
 * (c) Camille Scott, 2021
 * File   : diginorm.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 08.10.2021
 */

#ifndef GOETIA_DIGINORM_HH
#define GOETIA_DIGINORM_HH

#include "goetia/processors.hh"
#include "goetia/dbg.hh"
#include "goetia/storage/storage.hh"


namespace goetia {

template <class T>
struct DiginormFilter;

template <template <class, class> class GraphType,
                    class StorageType,
                    class ShifterType>
struct DiginormFilter<GraphType<StorageType, ShifterType>> {

    typedef GraphType<StorageType, ShifterType> graph_type;
    typedef ShifterType                         shifter_type;

    typedef typename shifter_type::alphabet     alphabet;
    typedef typename shifter_type::hash_type    hash_type;
   	typedef typename hash_type::value_type      value_type;
    typedef typename shifter_type::kmer_type    kmer_type;

    static bool median_count_at_least(const std::string& sequence,
                                      const unsigned int cutoff,
                                      graph_type&        graph) {

        auto kmers = graph.get_hash_iter(sequence);
        unsigned int min_req = 0.5 + float(sequence.size() - graph.K + 1) / 2;
        unsigned int num_cutoff_kmers = 0;

        // first loop:
        // accumulate at least min_req worth of counts before checking to see
        // if we have enough high-abundance k-mers to indicate success.
        for (unsigned int i = 0; i < min_req; ++i) {
            hash_type kmer = kmers->next();
            if (graph.query(kmer) >= cutoff) {
                ++num_cutoff_kmers;
            }
        }


        // second loop: now check to see if we pass the threshold for each k-mer.
        if (num_cutoff_kmers >= min_req) {
            return true;
        }
        while(!kmers->done()) {
            hash_type kmer = kmers->next();
            if (graph.query(kmer) >= cutoff) {
                ++num_cutoff_kmers;
                if (num_cutoff_kmers >= min_req) {
                    return true;
                }
            }
        }
        return false;
    }

    static auto median_in_range(const std::string& sequence,
                                const unsigned int low,
                                const unsigned int high,
                                graph_type&        graph)
    -> std::pair<bool, std::vector<hash_type>> {
        
        std::vector<goetia::count_t> counts;
        std::vector<hash_type> hashes;
        graph.query_sequence(sequence, counts, hashes);

        size_t n = counts.size() / 2;
        std::nth_element(counts.begin(), counts.begin() + n, counts.end());
        auto median = counts[n];
        bool in_range = (median > low && median <= high);

        return {in_range, hashes};
    }

    class Filter {
      private:
        
        std::shared_ptr<graph_type> graph;

      public:

        const uint16_t K;
        const unsigned int cutoff;

        Filter(std::shared_ptr<graph_type> graph,
               unsigned int                cutoff)
            : graph(graph),
              K(graph->K),
              cutoff(cutoff)
        {
        }

        static std::shared_ptr<Filter> build(std::shared_ptr<graph_type> graph,
                                             unsigned int                cutoff) {
            return std::make_shared<Filter>(graph, cutoff);
        }

        std::tuple<bool, uint64_t> filter_sequence(const std::string& sequence) {
            if (median_count_at_least(sequence, cutoff, *(graph.get()))) {
                return {false, sequence.length() - K + 1};
            }

            graph->insert_sequence(sequence);

            return {true, sequence.length() - K + 1};
        }

    };

    using Processor = FilterProcessor<Filter>;
};

}

#endif
