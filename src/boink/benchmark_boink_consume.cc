#include <iostream>

#include "boink/dbg.hh"
#include "boink/consumer.hh"

using namespace boink;

int main(int argc, char *argv[]) {
    DefaultDBG graph = DefaultDBG(21, 100000000, 4);
    FileConsumer<DefaultDBG> consumer = FileConsumer<DefaultDBG>(graph);
    std::string filename = std::string(argv[1]);

    uint64_t n_reads = 0;
    uint64_t n_kmers = 0;
    consumer.consume(filename, n_reads, n_kmers);

    std::cout << n_reads << " reads, " << n_kmers << " kmers consumed." << std::endl;
}
