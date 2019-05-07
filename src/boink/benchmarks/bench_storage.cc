#include "boink/benchmarks/bench_storage.hh"

#include "boink/storage/bitstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/nibblestorage.hh"
#include "boink/storage/sparseppstorage.hh"

#include <iostream>
#include <memory>
#include <random>

namespace boink {
namespace bench {



void run_storage_bench() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<hashing::hash_t> dis;

    std::vector<hashing::hash_t> hashes;
    for (int n = 0; n < 10000000; ++n) {
        hashes.push_back(dis(gen));
    }

    std::cout << "storage_type, bench, time" << std::endl;
    std::unique_ptr<storage::BitStorage> bitstorage = std::make_unique<storage::BitStorage>(5000000, 4);
    std::cout << "BitStorage, insert, " << time_it(storage_insert_bench<storage::BitStorage>, bitstorage, hashes) << std::endl
              << "BitStorage, query, " << time_it(storage_query_bench<storage::BitStorage>, bitstorage, hashes) << std::endl
              << "BitStorage, insert_and_query, " << time_it(storage_insert_and_query_bench<storage::BitStorage>, bitstorage, hashes) 
              << std::endl;

    std::unique_ptr<storage::NibbleStorage> nibblestorage = std::make_unique<storage::NibbleStorage>(5000000, 4);
    std::cout << "NibbleStorage, insert, " << time_it(storage_insert_bench<storage::NibbleStorage>, nibblestorage, hashes) 
              << ", query, " << time_it(storage_query_bench<storage::NibbleStorage>, nibblestorage, hashes) 
              << ", insert_and_query, " << time_it(storage_insert_and_query_bench<storage::NibbleStorage>, nibblestorage, hashes)
              << std::endl;

    std::unique_ptr<storage::ByteStorage> bytestorage = std::make_unique<storage::ByteStorage>(5000000, 4);
    std::cout << "ByteStorage, insert, " << time_it(storage_insert_bench<storage::ByteStorage>, bytestorage, hashes) << std::endl
              << "ByteStorage, query, " << time_it(storage_query_bench<storage::ByteStorage>, bytestorage, hashes) << std::endl 
              << "ByteStorage, insert_and_query, " << time_it(storage_insert_and_query_bench<storage::ByteStorage>, bytestorage, hashes) 
              << std::endl;

    std::unique_ptr<storage::SparseppSetStorage> sparseppstorage = std::make_unique<storage::SparseppSetStorage>();
    std::cout << "SparseppSetStorage, insert, " << time_it(storage_insert_bench<storage::SparseppSetStorage>, sparseppstorage, hashes) << std::endl
              << "SparseppSetStorage, query, " << time_it(storage_query_bench<storage::SparseppSetStorage>, sparseppstorage, hashes) << std::endl
              << "SparseppSetStorage, insert_and_query, " << time_it(storage_insert_and_query_bench<storage::SparseppSetStorage>, sparseppstorage, hashes) 
              << std::endl;

}

}
}
