/* cdbg.hh -- Streaming compact de Bruijn Graph
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

#define DEBUG_CDBG
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

typedef uint64_t id_t;
typedef uint64_t hash_t;
#define NULL_ID ULLONG_MAX
#define UNITIG_START_ID 100000000000

using std::make_shared;
using std::shared_ptr;
using std::string;
using namespace oxli;


typedef std::pair<hash_t, id_t> HashIDPair;
typedef std::unordered_set<hash_t> UHashSet;
typedef std::vector<hash_t> HashVector;
typedef std::unordered_map<hash_t, id_t> HashIDMap;
typedef std::unordered_set<id_t> IDSet;


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
    }
}



class CompactNode {
public:
    const id_t node_id;
    std::string sequence;
    
    CompactNode(id_t node_id, std::string sequence) :
        node_id(node_id), sequence(sequence) {}

    std::string revcomp() const {
        return _revcomp(sequence);
    }

    size_t length() const {
        return sequence.length();
    }

    friend bool operator== (const CompactNode& lhs, const CompactNode& rhs) {
        return lhs.node_id == rhs.node_id;
    }

};


class UnitigNode : public CompactNode {
public:

    id_t in_id, out_id;
    HashVector tags;


    UnitigNode(id_t node_id, std::string sequence,
               uint64_t length) :
        CompactNode(node_id, sequence),
        in_id(NULL_ID), out_id(NULL_ID) {}

};


class DecisionNode: public CompactNode {
public:
    using CompactNode::CompactNode;
    std::vector<id_t> in_edges;
    std::vector<id_t> out_edges;
    uint32_t count;

    DecisionNode(id_t node_id, std::string sequence) :
        CompactNode(node_id, sequence),
        count(0) {}

    uint8_t degree() const {
        return in_edges.size() + out_edges.size();
    }

    uint8_t in_degree() const {
        return in_edges.size();
    }

    uint8_t out_degree() const {
        return out_edges.size();
    }
};


typedef shared_ptr<CompactNode> CompactNodePtr;
typedef shared_ptr<DecisionNode> DecisionNodePtr;
typedef shared_ptr<UnitigNode> UnitigNodePtr;

class cDBG : public KmerClient {

    HashIDMap hash_id_map;
    std::unordered_map<id_t, CompactNodePtr> id_node_map;

public:

    cDBG(uint16_t K) :
        KmerClient(K) {

    }

};


template <class GraphType>
class StreamingCompactor : public AssemblerMixin<GraphType> {

    shared_ptr<GraphType> dbg;
    cDBG cdbg;

public:

    using AssemblerType = AssemblerMixin<GraphType>;
    using AssemblerType::seen;
    using AssemblerType::get_left;
    using AssemblerType::get_right;
    using AssemblerType::degree_left;
    using AssemblerType::degree_right;
    using AssemblerType::count_nodes;

    StreamingCompactor(shared_ptr<GraphType> dbg) :
        AssemblerMixin<GraphType>(dbg),
        dbg(dbg),
        cdbg(dbg->K()) {

    }

    std::string compactify(const std::string seed) {
        Path path;
        this->set_cursor(seed);
        this->get_cursor(path);
        compactify_left(path);
        this->set_cursor(seed);
        compactify_right(path);

        return this->to_string(path);
    }

    void compactify_right(Path& path) {
        this->seen.insert(this->get());
        
        shift_t next;
        while (get_right(next)
               && !this->seen.count(next.hash)
               && !(degree_left() > 1)) {
            path.push_back(next.symbol);
            this->seen.insert(next.hash);
        }
    }

    void compactify_left(Path& path) {
        this->seen.insert(this->get());
        
        shift_t next;
        while (get_left(next)
               && !this->seen.count(next.hash)
               && !(degree_right() > 1)) {
            path.push_front(next.symbol);
            this->seen.insert(next.hash);
        }

    }

    bool is_decision_node(uint8_t& degree) {
        // oops better put notes here
        uint8_t ldegree, rdegree;
        ldegree = this->degree_left();
        rdegree = this->degree_right();
        degree = ldegree + rdegree;
        return ldegree > 1 || rdegree > 1;
    }

    bool is_decision_node(const std::string& node,
                          uint8_t& degree) {
        this->set_cursor(node);
        return is_decision_node(degree);
    }

    bool is_decision_node(const std::string& node) {
        this->set_cursor(node);
        return this->degree_left() > 1 || this->degree_right() > 1;
    }

    bool insert_sequence(const std::string& sequence,
                         std::vector<uint32_t>& decision_positions,
                         HashVector& decision_hashes,
                         std::vector<NeighborBundle>& decision_neighbors) {

        if (!dbg->add_sequence(sequence)) {
            return false; // if there were no new k-mers, nothing to do
        }

        KmerIterator<typename GraphType::shifter_type> iter(sequence, this->_K);
        size_t pos = 0;
        while(!iter.done()) {
            hash_t h = iter.next();
            std::vector<shift_t> left_neighbors = iter.shifter.gather_left();
            std::vector<shift_t> right_neighbors = iter.shifter.gather_right();

            if (count_nodes(left_neighbors) > 1 || count_nodes(right_neighbors) > 1) {
                StringVector left_kmers, right_kmers;
                for (auto neighbor : left_neighbors) {
                    left_kmers.push_back(neighbor.symbol
                                         + sequence.substr(pos, this->_K-1));
                }
                for (auto neighbor : right_neighbors) {
                    right_kmers.push_back(sequence.substr(pos + 1, this->_K-1)
                                          + neighbor.symbol);
                                          
                }
                decision_positions.push_back(pos);
                decision_hashes.push_back(h);
                decision_neighbors.push_back(std::make_pair(left_kmers, right_kmers));
            }

            ++pos;
        }

        return true;
    }

    /*
    void update(const std::string& sequence) {
        std::vector<uint32_t> decision_positions;
        std::vector<hash_t> decision_hashes;
        std::vector<NeighborBundle> decision_neighbors;
        if (insert_sequence(sequence,
                            decision_positions,
                            decision_hashes,
                            decision_neighbors)) {

            
            if (decision_positions.size() == 0) {
                // sequence doesn't have an decision nodes,
                // do a "base case" upate
                update_linear(sequence);
            } else {
                
            }

        } else {
            return;
        }
    }

    void update(const std::string& sequence) {
        if (!dbg->get(seed)) {
            return;
        } else if (is_decision_node(seed)) {
            // find neighbors and update from them
        } else {
            // update unitig
            Path unitig;
            this->set_cursor(seed);
            
        }
    }
    */

};



}


#endif
