#include <iostream>

#include "oxli/hashtable.hh"

using namespace oxli;
using namespace oxli::read_parsers;

int main(int argc, char *argv[]) {
    std::vector<uint64_t> sizes = get_n_primes_near_x(4, 100000000);
    Nodetable graph = Nodetable(21, sizes);
    std::string filename = std::string(argv[1]);

    unsigned int n_reads = 0;
    long long unsigned int  n_kmers = 0;
    graph.consume_seqfile<FastxReader>(filename, n_reads, n_kmers);

    std::cout << n_reads << " reads, " << n_kmers << " kmers consumed." << std::endl;
}
