#include <iostream>

#include "goetia/dbg.hh"
#include "goetia/storage/bitstorage.hh"
#include "goetia/hashing/rollinghashshifter.hh"
#include "goetia/processors.hh"

using namespace goetia;

typedef dBG<BitStorage, hashing::DefaultShifter> GraphType;

int main(int argc, char *argv[]) {
    auto graph = std::make_shared<GraphType>(21, 100000000, 4);
    FileConsumer<GraphType> consumer(graph);
    std::string filename = std::string(argv[1]);
    
    consumer.process(filename);
}
