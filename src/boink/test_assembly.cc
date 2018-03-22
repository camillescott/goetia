#include "boink/assembly.hh"
#include "boink/boink.hh"
#include "boink/dbg.hh"

#include <cassert>
#include <string>
#include <iostream>

using namespace boink;
namespace boink {
void test_assemble_right() {
    uint16_t K = 21;
    auto g = make_shared<DefaultDBG>(21, 1000, 4);
    auto assembler = make_assembler(g);

    string S = "TACATGCTTACTCACAGCACCCCTTCTAATGCGTAACCAGGCGTGGAATT";
    g->add_sequence(S);
    std::cout << "added " << S << std::endl;
    auto counts = g->get_counts(S);
    for (auto c : counts) {
        std::cout << std::to_string(c) << " ";
    }
    std::cout << std::endl;
    
    Path path;

    assembler.assemble(S.substr(0, K), path);
    auto res = std::string(path.begin(), path.end());
    std::cout << "assemble from " << S.substr(0, K) << std::endl;
    std::cout << "expected:  " << S << std::endl;
    std::cout << "assembled: " << res << std::endl;
    
    assert(S == res);
}

void test_assemble_left() {
    uint16_t K = 21;
    auto g = make_shared<DefaultDBG>(21, 1000, 4);
    auto assembler = make_assembler(g);

    string S = "TACATGCTTACTCACAGCACCCCTTCTAATGCGTAACCAGGCGTGGAATT";
    g->add_sequence(S);

    Path path;
    assembler.assemble(S.substr(S.length()-K), path);
    auto res = std::string(path.begin(), path.end());
    std::cout << "assemble from " << S.substr(S.length()-K) << std::endl;
    std::cout << "expected:  " << S << std::endl;
    std::cout << "assembled: " << res << std::endl;

    assert(S == res);
}

}

int main() {
    boink::test_assemble_right();
    boink::test_assemble_right();
}
