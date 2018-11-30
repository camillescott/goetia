#include "boink/hashing/alphabets.hh"
#include <iostream>

using namespace boink;

int main(int argc, char *argv[]) {
    std::string kmer = "AAAAAAAAAA";
    std::cout << hashing::hash_cyclic(kmer, 10) << std::endl;
}
