#include "goetia/benchmarks/bench_storage.hh"

#include "goetia/storage/bitstorage.hh"
#include "goetia/storage/bytestorage.hh"
#include "goetia/storage/nibblestorage.hh"
#include "goetia/storage/sparseppstorage.hh"
#include "goetia/storage/phmapstorage.hh"

#include <iostream>
#include <memory>
#include <random>

namespace goetia {
namespace bench {

std::vector<uint64_t> generate_hashes(size_t n_hashes) {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis;

    std::vector<uint64_t> hashes;
    for (int n = 0; n < n_hashes; ++n) {
        hashes.push_back(dis(gen));
    }

    return hashes;
}


void run_storage_bench() {

    std::vector<size_t> hashes_sizes = {1000000, 10000000, 100000000};


    std::unique_ptr<storage::BitStorage> bitstorage;
    std::unique_ptr<storage::NibbleStorage> nibblestorage;
    std::unique_ptr<storage::ByteStorage> bytestorage;
    std::unique_ptr<storage::SparseppSetStorage> sparseppstorage;
    std::unique_ptr<storage::PHMapStorage> phmapstorage;
    
    std::cout << "storage_type, n_hashes, bench, time" << std::endl;
    for (auto n_hashes : hashes_sizes) {
        bitstorage = std::make_unique<storage::BitStorage>(n_hashes / 4, 4);
        nibblestorage = std::make_unique<storage::NibbleStorage>(n_hashes / 4, 4);
        bytestorage = std::make_unique<storage::ByteStorage>(n_hashes / 4, 4);
        sparseppstorage  = std::make_unique<storage::SparseppSetStorage>();
        phmapstorage = std::make_unique<storage::PHMapStorage>();

        auto hashes = generate_hashes(n_hashes);

        for (size_t N = 0; N < 3; ++N) {
            _run_storage_bench(bitstorage, hashes, "BitStorage");
            _run_storage_bench(nibblestorage, hashes, "NibbleStorage");
            _run_storage_bench(bytestorage, hashes, "ByteStorage");
            _run_storage_bench(phmapstorage, hashes, "PHMapStorage");
            _run_storage_bench(sparseppstorage, hashes, "SparseppSetStorage");
        }
    }
}

}
}
