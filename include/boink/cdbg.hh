/* cdbg.hh -- compact de Bruijn Graph
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef CDBG_HH
#define CDBG_HH

#include <algorithm>
#include <cstdint>
#include <functional>
#include <memory>
#include <limits>
#include <list>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>

#include "oxli/alphabets.hh"
#include "boink/assembly.hh"
#include "boink/hashing.hh"
#include "boink/dbg.hh"
#include "boink/minimizers.hh"

# ifdef DEBUG_CDBG
#   define pdebug(x) do { std::cerr << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif

#define complement(ch) ((ch) == 'A' ? 'T' : \
                        (ch) == 'T' ? 'A' : \
                        (ch) == 'C' ? 'G' : 'C')

namespace boink {

using std::string;
using std::unique_ptr;
using std::make_unique;
using std::vector;
using std::pair;

using namespace oxli;

typedef uint64_t id_t;
typedef pair<hash_t, hash_t> junction_t;


// hash_combine and pair_hash courtesy SO:
// https://stackoverflow.com/questions/9729390/how-to-use-
// unordered-set-with-custom-types/9729747#9729747
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
}


struct junction_hash
{
    inline std::size_t operator()(const junction_t & v) const
    {
         std::size_t seed = 0;
         hash_combine(seed, v.first);
         hash_combine(seed, v.second);
         return seed;
    }
};

template<typename _Ty1, typename _Ty2>
std::ostream& operator<<(std::ostream& _os, const std::pair<_Ty1, _Ty2>& _p) {
    _os << "(" << _p.first << ", " << _p.second << ")";
    return _os;
}


#define NULL_ID             ULLONG_MAX
#define UNITIG_START_ID     0
#define NULL_JUNCTION       make_pair(0,0)

typedef vector<hash_t> HashVector;


enum node_meta_t {
    FULL,
    TIP,
    ISLAND,
    TRIVIAL
};


inline const char * node_meta_repr(node_meta_t meta) {
    switch(meta) {
        case FULL:
            return "FULL";
        case TIP:
            return "TIP";
        case ISLAND:
            return "ISLAND";
        case TRIVIAL:
            return "TRIVIAL";
        default:
            return "UNKNOWN";
    }
}


struct node_meta_counter {

    int64_t full_count;
    int64_t tip_count;
    int64_t island_count;
    int64_t unknown_count;
    int64_t trivial_count;

    node_meta_counter() :
        full_count(0),
        tip_count(0),
        island_count(0),
        unknown_count(0),
        trivial_count(0) {
    }

    void mutate(node_meta_t meta, int64_t amt) {
        switch(meta) {
            case FULL:
                full_count += amt;
                break;
            case TIP:
                tip_count += amt;
                break;
            case ISLAND:
                island_count += amt;
                break;
            case TRIVIAL:
                trivial_count += amt;
                break;
            default:
                unknown_count += amt;
                break;
        }
    }

    void increment(node_meta_t meta) {
        mutate(meta, 1);
    }

    void decrement(node_meta_t meta) {
        mutate(meta, -1);
    }

    string repr() const {
        std::ostringstream os;
        os << *this;
        return os.str();
    }

    string header() const {
        return string("full,tip,island,trivial,unknown");
    }

    friend std::ostream& operator<<(std::ostream& o, const node_meta_counter& c);

};

std::ostream& operator<<(std::ostream& o, const node_meta_counter& c) {
    o << c.full_count << "," << c.tip_count << "," << c.island_count
      << "," << c.trivial_count << "," << c.unknown_count;
    return o;
}


class CompactNode {
public:
    const id_t node_id;
    string sequence;
    
    CompactNode(id_t node_id, const string& sequence) :
        node_id(node_id), sequence(sequence) {}

    string revcomp() const {
        return _revcomp(sequence);
    }

    size_t length() const {
        return sequence.length();
    }

    friend bool operator== (const CompactNode& lhs, const CompactNode& rhs) {
        return lhs.node_id == rhs.node_id;
    }

};


class DecisionNode: public CompactNode {

protected:

    bool _dirty;

public:
    vector<junction_t> left_juncs;
    vector<junction_t> right_juncs;
    uint32_t count;

    DecisionNode(id_t node_id, const string& sequence) :
        CompactNode(node_id, sequence),
        _dirty(true),
        count(1) {
        
    }

    bool is_dirty() const {
        return _dirty;
    }

    void set_dirty(bool dirty) {
        _dirty = dirty;
    }

    const size_t degree() const {
        return left_degree() + right_degree();
    }

    const size_t left_degree() const {
        return left_juncs.size();
    }

    const size_t right_degree() const {
        return right_juncs.size();
    }

    bool has_left_junc(junction_t j) const {
        for (auto junc : left_juncs) {
            if (j == junc) return true;
        }
        return false;
    }

    void add_left_junc(junction_t j) {
        if (!has_left_junc(j)) {
            left_juncs.push_back(j);
        }
    }

    void remove_left_junc(junction_t j) {
        for (auto it = left_juncs.begin(); it != left_juncs.end(); ) {
            if (*it == j) {
                left_juncs.erase(it);
                return;
            }
        }
    }

    bool has_right_junc(junction_t j) const {
        for (auto junc : right_juncs) {
            if (j == junc) return true;
        }
        return false;
    }

    void add_right_junc(junction_t j) {
        if (!has_right_junc(j)) {
            right_juncs.push_back(j);
        }
    }

    void remove_right_junc(junction_t j) {
        for (auto it = right_juncs.begin(); it != right_juncs.end(); ) {
            if (*it == j) {
                right_juncs.erase(it);
                return;
            }
        }
    }

    std::string repr() const {
        std::ostringstream os;
        os << *this;
        return os.str();
    }

    friend std::ostream& operator<<(std::ostream& o, const DecisionNode& dn);
};


std::ostream& operator<<(std::ostream& o, const DecisionNode& dn) {

    o << "<DNode ID/hash=" << dn.node_id << " k-mer=" << dn.sequence
      << " Dl=" << std::to_string(dn.left_degree())
      << " Dr=" << std::to_string(dn.right_degree())
      << " count=" << dn.count << " dirty=" << dn.is_dirty() << ">";
    return o;
}


class UnitigNode : public CompactNode {
public:

    junction_t left_junc, right_junc;
    HashVector tags;
    node_meta_t meta;

    UnitigNode(id_t node_id,
               junction_t left_junc,
               junction_t right_junc,
               const string& sequence)
        : CompactNode(node_id, sequence),
          left_junc(left_junc),
          right_junc(right_junc) {
        
    }

    std::string repr() const {
        std::ostringstream os;
        os << *this;
        return os.str();
    }

    friend std::ostream& operator<<(std::ostream& o, const UnitigNode& un);

};


std::ostream& operator<<(std::ostream& o, const UnitigNode& un) {
    o << "<UNode ID=" << un.node_id
      << " left=(" << un.left_junc.first 
      << "," << un.left_junc.second << ")"
      << " right=(" << un.right_junc.first 
      << "," << un.right_junc.second << ")"
      << " meta=" << node_meta_repr(un.meta)
      << ">";
    return o;
}


typedef CompactNode * CompactNodePtr;
typedef DecisionNode * DecisionNodePtr;
typedef UnitigNode * UnitigNodePtr;

class cDBG : public KmerClient {

public:

    typedef std::unordered_map<hash_t,
                               std::unique_ptr<DecisionNode>> dnode_map_t;
    typedef dnode_map_t::const_iterator dnode_iter_t;
    typedef std::unordered_map<id_t,
                               std::unique_ptr<UnitigNode>> unode_map_t;
    typedef unode_map_t::const_iterator unode_iter_t;

protected:

    dnode_map_t decision_nodes;
    unode_map_t unitig_nodes;
    std::unordered_map<junction_t, id_t, junction_hash> unitig_junction_map;
    std::unordered_map<hash_t, id_t> unitig_tag_map;

    uint64_t _n_updates;
    uint64_t _unitig_id_counter;
    uint64_t _n_unitig_nodes;

public:

    enum Graph_Ops_t {
        ADD_DNODE,
        ADD_UNODE,
        DELETE_DNODE,
        DELETE_UNODE,
        ADD_EDGE,
        INCR_DNODE_COUNT
    };

    node_meta_counter meta_counter;

    cDBG(uint16_t K) :
        KmerClient(K),
        _n_updates(0),
        _unitig_id_counter(UNITIG_START_ID),
        _n_unitig_nodes(0) {

    }

    dnode_map_t::const_iterator dnodes_begin() const {
        return decision_nodes.cbegin();
    }

    dnode_map_t::const_iterator dnodes_end() const {
        return decision_nodes.cend();
    }

    unode_map_t::const_iterator unodes_begin() const {
        return unitig_nodes.cbegin();
    }

    unode_map_t::const_iterator unodes_end() const {
        return unitig_nodes.cend();
    }

    uint64_t n_updates() const {
        return _n_updates;
    }

    uint64_t n_unitig_nodes() const {
        return _n_unitig_nodes;
    }

    uint64_t n_decision_nodes() const {
        return decision_nodes.size();
    }

    uint64_t n_tags() const {
        return unitig_tag_map.size();
    }

    DecisionNode* build_dnode(hash_t hash, const string& kmer) {
        DecisionNode * dnode = get_dnode(hash);
        if (dnode == nullptr) {
            pdebug("Build d-node " << hash << ", " << kmer);
            unique_ptr<DecisionNode> dnode_ptr = make_unique<DecisionNode>(hash, kmer);
            decision_nodes.insert(make_pair(hash, std::move(dnode_ptr)));
            dnode = get_dnode(hash);
        }
        return dnode;
    }

    DecisionNode* get_dnode(hash_t hash) {
        auto search = decision_nodes.find(hash);
        if (search != decision_nodes.end()) {
            return search->second.get();
        }
        return nullptr;
    }

    template<class ShifterType>
    vector<DecisionNode*> get_dnodes(const string& sequence) {
        KmerIterator<ShifterType> iter(sequence, this->_K);
        vector<DecisionNode*> result;
        while(!iter.done()) {
            hash_t h = iter.next();
            DecisionNode * dnode;
            if ((dnode = get_dnode(h)) != nullptr) {
                result.push_back(dnode);
            }
        }
        
        return result;
    }

    UnitigNode * build_unode(HashVector& tags,
                             const string& sequence,
                             junction_t left_junc,
                             junction_t right_junc) {

        pdebug("Attempt u-node build on..."
                << " left=(" << left_junc.first 
                << "," << left_junc.second << ")"
                << " right=(" << right_junc.first 
                << "," << right_junc.second << ")");
        // Check for existing Unitigs on these junctions
        UnitigNode * existing_left, * existing_right;
        bool left_valid = false, right_valid = false;
        existing_left = get_unode(left_junc);
        if (existing_left) {
            pdebug("Found existing from left_junc " << *existing_left);
            if (left_junc == existing_left->left_junc &&
                right_junc == existing_left->right_junc) {
                left_valid = true;
            } else {
                delete_unode(existing_left);
            }
        }
        existing_right = get_unode(right_junc);
        if (existing_right) {
            pdebug("Found existing from right_junc " << *existing_right);
            if (left_junc == existing_right->left_junc &&
                right_junc == existing_right->right_junc) {
                right_valid = true;
            } else {
                delete_unode(existing_right);
            }
        }

        if (left_valid) {
            return existing_left;
        }
        if (right_valid) {
            return existing_right;
        }

        // No valid existing Unitigs, make a new one
        id_t id = _unitig_id_counter;
        unique_ptr<UnitigNode> unode = make_unique<UnitigNode>(id,
                                                               left_junc,
                                                               right_junc,
                                                               sequence);
        _unitig_id_counter++;
        _n_unitig_nodes++;

        // Link up its new tags
        unode->tags.insert(std::end(unode->tags), std::begin(tags), std::end(tags));
        for (auto tag: tags) {
            unitig_tag_map.insert(make_pair(tag, id));
        }

        // Link up its junctions
        unitig_junction_map.insert(make_pair(left_junc, id));
        unitig_junction_map.insert(make_pair(right_junc, id));
        link_unode_to_dnodes(unode.get());

        unode->meta = get_unode_meta(unode.get());
        meta_counter.increment(unode->meta);
        pdebug("Built unode " << *(unode.get()));

        // Transfer the UnitigNode's ownership to the map;
        // get its new memory address and return it
        unitig_nodes.insert(make_pair(id, std::move(unode)));
        return unitig_nodes[id].get();
    }

    void link_unode_to_dnodes(UnitigNode * unode) {
        DecisionNode * left = get_left_dnode(unode);
        DecisionNode * right = get_right_dnode(unode);
        if (left != nullptr) {
            left->add_right_junc(unode->left_junc);
        }
        if (right != nullptr) {
            right->add_left_junc(unode->right_junc);
        }
    }

    DecisionNode * get_left_dnode(UnitigNode * unode) {
        return get_dnode(unode->left_junc.first);
    }

    DecisionNode * get_right_dnode(UnitigNode * unode) {
        return get_dnode(unode->right_junc.second);
    }

    node_meta_t get_unode_meta(UnitigNode * unode) {
        bool has_left = false, has_right = false;
        has_left = (get_left_dnode(unode) != nullptr);
        has_right = (get_right_dnode(unode) != nullptr);
        if (has_left && has_right) {
            if (unode->left_junc == unode->right_junc) {
                return TRIVIAL;
            } else {
                return FULL;
            }
        } else if (has_left != has_right) {
            return TIP;
        } else {
            return ISLAND;
        }
    }

    UnitigNode * get_unode(junction_t junc) {
        auto search = unitig_junction_map.find(junc);
        if (search != unitig_junction_map.end()) {
            id_t id = search->second;
            return unitig_nodes[id].get();
        }
        return nullptr;
    }

    UnitigNode * get_unode(hash_t hash) {
        auto search = unitig_tag_map.find(hash);
        if (search != unitig_tag_map.end()) {
            id_t id = search->second;
            return unitig_nodes[id].get();
        }
        return nullptr;
    }

    UnitigNode * get_unode_from_id(id_t id) {
        auto search = unitig_nodes.find(id);
        if (search != unitig_nodes.end()) {
            return search->second.get();
        }
        return nullptr;
    }

    void delete_unode(UnitigNode * unode) {
        if (unode != nullptr) {
            pdebug("Deleting " << *unode);
            id_t id = unode->node_id;
            meta_counter.decrement(unode->meta);
            for (hash_t tag: unode->tags) {
                unitig_tag_map.erase(tag);
            }
            unitig_junction_map.erase(unode->left_junc);
            unitig_junction_map.erase(unode->right_junc);

            unitig_nodes.erase(id);
            unode = nullptr;
            _n_unitig_nodes--;
            _n_updates++;
        }
    }

    void write_adj_matrix(const string& filename) {
        std::ofstream out;
        out.open(filename);

        std::cerr << "Gather node IDs." << std::endl;
        std::vector<id_t> dnode_ids(decision_nodes.size());
        std::vector<id_t> unode_ids(unitig_nodes.size());
        size_t i = 0;
        for (auto it = decision_nodes.begin(); it != decision_nodes.end(); ++it) {
            dnode_ids[i] = it->first;
            ++i;
        }
        i = 0;
        for (auto it = unitig_nodes.begin(); it != unitig_nodes.end(); ++it) {
            unode_ids[i] = it->first;
            ++i;
        }
        std::cerr << "Sorting nodes." << std::endl;
        std::sort(dnode_ids.begin(), dnode_ids.end());
        std::sort(unode_ids.begin(), unode_ids.end());

        std::cerr << "Writing out..." << std::endl;

        // one row at a time
        for (i = 0; i < dnode_ids.size() + unode_ids.size(); ++i) {
            std::set<id_t> neighbors;
            if (i < dnode_ids.size()) {
                id_t root_id = dnode_ids[i];
                for (auto junc : decision_nodes[root_id]->left_juncs) {
                    id_t neighbor_id = unitig_junction_map[junc];
                    neighbors.insert(neighbor_id);
                }
                for (auto junc : decision_nodes[root_id]->right_juncs) {
                    id_t neighbor_id = unitig_junction_map[junc];
                    neighbors.insert(neighbor_id);
                }
                for (size_t j = 0; j < dnode_ids.size(); ++j) {
                    out << "0 ";
                }
                for (auto neighbor_id : unode_ids) {
                    if (neighbors.count(neighbor_id)) {
                        out << neighbor_id;
                    } else {
                        out << "0";
                    }
                    out << " ";
                }
            } else {
                id_t root_id = unode_ids[i];
                id_t left = unitig_nodes[root_id]->left_junc.first;
                id_t right = unitig_nodes[root_id]->right_junc.second;
                if (decision_nodes.count(left)) {
                    neighbors.insert(left);
                }
                if (decision_nodes.count(right)) {
                    neighbors.insert(right);
                }
                for (auto neighbor_id : dnode_ids) {
                    if (neighbors.count(neighbor_id)) {
                        out << neighbor_id;
                    } else {
                        out << "0";
                    }
                    out << " ";
                }
                for (size_t j = 0; j < unode_ids.size(); ++j) {
                    out << "0 ";
                }
            }

            out << std::endl;
        }

        std::cerr << "Wrote " << i << "x" << i << " adjacency matrix." << std::endl;
    }
};
}

#undef pdebug
#endif
