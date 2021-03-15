/* (c) Camille Scott, 2021
 * File   : phmapstorage.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 12.03.2021
 */

#include "goetia/goetia.hh"
#include "goetia/storage/btreestorage.hh"

#include <cstdint>
#include <cstring>
#include <sstream> // IWYU pragma: keep

namespace goetia {

namespace storage {


const count_t
BTreeStorage::insert_and_query(value_type h) {
    insert(h);
    return 1; // its a presence filter so always 1 after insert
}


const count_t
BTreeStorage::query(value_type h) const {
    return _store->count(h);
}


std::shared_ptr<BTreeStorage>
BTreeStorage::build() {
    return std::make_shared<BTreeStorage>();
}


std::shared_ptr<BTreeStorage>
BTreeStorage::build(const typename StorageTraits<BTreeStorage>::params_type&) {
    return std::make_shared<BTreeStorage>();
}


std::shared_ptr<BTreeStorage>
BTreeStorage::clone() const {
    return std::make_shared<BTreeStorage>();
}


void BTreeStorage::save(std::string filename, uint16_t K) {

}

void BTreeStorage::load(std::string filename, uint16_t &K) {

}

void BTreeStorage::serialize(std::ofstream& out) {
    //out.write(std::string(this->NAME).c_str(), this->NAME.size());
    //out.write(this->version_binary(), sizeof(this->OBJECT_ABI_VERSION));
    //_store->serialize(BaseSppSerializer(), &out);
}

std::shared_ptr<BTreeStorage>
BTreeStorage::deserialize(std::ifstream& in) {
    /*
    std::string name;
    name.resize(Tagged<BTreeStorage>::NAME.size());
    size_t version;

    in.read(name.data(), name.size());
    in.read(reinterpret_cast<char *>(&version), sizeof(version));

    if (name != Tagged<BTreeStorage>::NAME) {
        std::ostringstream err;
        err << "File has wrong type tag: found "
            << name
            << ", should be "
            << Tagged<BTreeStorage>::NAME;
        throw GoetiaFileException(err.str());
    } else if (version != Tagged<BTreeStorage>::OBJECT_ABI_VERSION) {
        std::ostringstream err;
        err << "File has wrong binary version: found "
            << std::to_string(version)
            << ", expected "
            << std::to_string(Tagged<BTreeStorage>::OBJECT_ABI_VERSION);
        throw GoetiaFileException(err.str());

    }

    auto storage = BTreeStorage::build();
    storage->_store->unserialize(BaseSppSerializer(), &in);
    return storage;
    */
}

}

}
