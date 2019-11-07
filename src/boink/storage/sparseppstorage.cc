/**
 * (c) Camille Scott, 2019
 * File   : sparseppstorage.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 30.08.2019
 */

#include "boink/boink.hh"
#include "boink/storage/sparseppstorage.hh"
#include "boink/storage/sparsepp/spp.h"
#include "boink/storage/sparsepp/serialize.hh"

#include <cstdint>
#include <sstream> // IWYU pragma: keep

namespace boink {

template <> const std::string Tagged<storage::SparseppSetStorage>::NAME = "Boink::SparseppSetStorage";

namespace storage {

const bool
SparseppSetStorage::insert(value_type h) {
    auto result = _store->insert(h);
    // the second in the returned pair reports that the insert
    // took place ie the hash was new
    return result.second;
}


const count_t
SparseppSetStorage::insert_and_query(value_type h) {
    insert(h);
    return 1; // its a presence filter so always 1 after insert
}


const count_t
SparseppSetStorage::query(value_type h) const {
    return _store->count(h);
}


std::shared_ptr<SparseppSetStorage>
SparseppSetStorage::build() {
    return std::make_shared<SparseppSetStorage>();
}


std::shared_ptr<SparseppSetStorage>
SparseppSetStorage::clone() const {
    return std::make_shared<SparseppSetStorage>();
}


void SparseppSetStorage::save(std::string filename, uint16_t K) {

}

void SparseppSetStorage::load(std::string filename, uint16_t &K) {

}

void SparseppSetStorage::serialize(std::ofstream& out) {
    out.write(this->NAME.c_str(), this->NAME.size());
    out.write(reinterpret_cast<const char*>(&(this->VERSION)), sizeof(this->VERSION));
    _store->serialize(BaseSppSerializer(), &out);
}

std::shared_ptr<SparseppSetStorage> SparseppSetStorage::deserialize(std::ifstream& in) {
    std::string name;
    name.resize(Tagged<SparseppSetStorage>::NAME.size());
    unsigned int version;

    in.read(name.data(), name.size());
    in.read(reinterpret_cast<char *>(&version), sizeof(version));

    if (name != Tagged<SparseppSetStorage>::NAME) {
        std::ostringstream err;
        err << "File has wrong type tag: found "
            << name
            << ", should be "
            << Tagged<SparseppSetStorage>::NAME;
        throw BoinkFileException(err.str());
    } else if (version != Tagged<SparseppSetStorage>::VERSION) {
        std::ostringstream err;
        err << "File has wrong binary version: found "
            << std::to_string(version)
            << ", expected "
            << std::to_string(Tagged<SparseppSetStorage>::VERSION);
        throw BoinkFileException(err.str());

    }

    //Tagged<SparseppSetStorage>::VERSION;
    auto storage = SparseppSetStorage::build();
    storage->_store->unserialize(BaseSppSerializer(), &in);
    return storage;
}

}

}
