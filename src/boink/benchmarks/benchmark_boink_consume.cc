#include <iostream>

#include "boink/dbg.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/processors.hh"

using namespace boink;

typedef dBG<storage::BitStorage, hashing::DefaultShifter> GraphType;

int main(int argc, char *argv[]) {
    auto graph = std::make_shared<GraphType>(21, 100000000, 4);
    FileConsumer<GraphType> consumer(graph);
    std::string filename = std::string(argv[1]);
    
    consumer.process(filename);
}
