#include "goetia/hashing/hashextender.hh"

#include <iostream>

using namespace goetia;

int main() {
    std::string known_kmer = "TCACCTGTGTTGTGCTACTTGCGGCGC";
    uint64_t known_hash =  13194817695400542713;

    FwdRollingExtender from_k(27);
    FwdRollingExtender from_other(from_k);

    std::cout << "=== from k ===" << std::endl;
    std::cout << from_k.set_cursor(known_kmer) << std::endl << from_k.get() << std::endl;

    std::cout << "=== from other ===" << std::endl;
    std::cout << from_other.set_cursor(known_kmer) << std::endl << from_other.get() << std::endl;

    CanRollingExtender can_from_k(27);
    CanRollingExtender can_from_other(can_from_k);

    std::cout << "=== from k ===" << std::endl;
    std::cout << can_from_k.set_cursor(known_kmer) << std::endl << can_from_k.get() << std::endl;

    std::cout << "=== from other ===" << std::endl;
    std::cout << can_from_other.set_cursor(known_kmer) << std::endl << can_from_other.get() << std::endl;
}
