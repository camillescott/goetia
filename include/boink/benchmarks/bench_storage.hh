#include "boink/boink.hh"
#include "boink/storage/storage.hh"
#include "boink/hashing/hashing_types.hh"

#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <vector>


namespace boink {
namespace bench {

template<class storage_type>
void storage_insert_bench(std::unique_ptr<storage_type>& storage,
                          std::vector<hashing::hash_t>& hashes) {

    for (auto hash: hashes) {
        storage->insert(hash);
    }
}

template<class storage_type>
void storage_query_bench(std::unique_ptr<storage_type>& storage,
                         std::vector<hashing::hash_t>& hashes) {

    for (auto hash: hashes) {
        auto val = storage->query(hash);
    }
}




template<class storage_type>
void storage_insert_and_query_bench(std::unique_ptr<storage_type>& storage,
                                    std::vector<hashing::hash_t>& hashes)  {

    for (auto hash: hashes) {
        auto val = storage->insert_and_query(hash);
    }
}


template<class callable, class storage_type>
double time_it(callable &&func,
               std::unique_ptr<storage_type>& storage,
               std::vector<hashing::hash_t>& hashes)  {
    
    auto time_start = std::chrono::system_clock::now();
    func(storage, hashes);
    auto time_elapsed = std::chrono::system_clock::now() - time_start;
    return std::chrono::duration<double>(time_elapsed).count();
}


template<class storage_type>
void _run_storage_bench(std::unique_ptr<storage_type>& storage,
                   std::vector<hashing::hash_t>& hashes,
                   std::string storage_name) {

    storage->reset();
    std::cout << storage_name << ", " << hashes.size() << ", insert, "           << time_it(storage_insert_bench<storage_type>, storage, hashes) << std::endl;
    std::cout << storage_name << ", " << hashes.size() << ", insert_second, "    << time_it(storage_insert_bench<storage_type>, storage, hashes) << std::endl;
    std::cout << storage_name << ", " << hashes.size() << ", query, "            << time_it(storage_query_bench<storage_type>, storage, hashes) << std::endl;
    std::cout << storage_name << ", " << hashes.size() << ", insert_and_query, " << time_it(storage_insert_and_query_bench<storage_type>, storage, hashes) << std::endl;
}

std::vector<hashing::hash_t> generate_hashes(size_t n_hashes);


void run_storage_bench();

}
}
