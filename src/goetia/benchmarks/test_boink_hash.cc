#include "goetia/hashing/alphabets.hh"
#include <iostream>

using namespace goetia;

int main(int argc, char *argv[]) {
    std::string kmer = "AAAAAAAAAA";
    std::cout << hashing::hash_cyclic(kmer, 10) << std::endl;
}
