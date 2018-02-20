
#include "oxli/hashtable.hh"
#include "boink/dbg.hh"

using namespace boink;

namespace boink {

template <class StorageType,
          class HashShifter>
std::vector<bool> dBG<StorageType, HashShifter>::
    add_sequence(const std::string& sequence) {

    KmerIterator<HashShifter> iter(sequence, _K);
    std::vector<bool> consumed(sequence.length() - _K + 1);

    size_t pos = 0;
    while(!iter.done()) {
        hash_t h = iter.next();
        consumed[pos] = add(h);
        ++pos;
    }

    return consumed;
}

};

int main() {
    std::string S = "ATCGGGCTGGGGTGCGCGGCATTATAAGGCTGGATG";
    boink::DefaultDBG G(16, oxli::get_n_primes_near_x(4, 100));

    G.add_sequence(S);
}
