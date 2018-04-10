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

using std::string;
using std::unique_ptr;
using std::make_unique;
using std::vector;
using std::pair;

using namespace oxli;


typedef pair<hash_t, id_t> HashIDPair;
typedef std::unordered_set<hash_t> UHashSet;
typedef vector<hash_t> HashVector;
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


class UnitigNode : public CompactNode {
public:

    id_t in_id, out_id;
    HashVector tags;

    UnitigNode(id_t node_id, const string& sequence) :
        CompactNode(node_id, sequence),
        in_id(NULL_ID), out_id(NULL_ID) {}

};


class DecisionNode: public CompactNode {

protected:

    bool _dirty;

public:
    using CompactNode::CompactNode;
    vector<id_t> in_edges;
    vector<id_t> out_edges;
    uint32_t count;

    DecisionNode(id_t node_id, const string& sequence) :
        CompactNode(node_id, sequence),
        count(1),
        _dirty(true) {
        
    }

    bool is_dirty() const {
        return _dirty;
    }

    uint8_t degree() const {
        return in_edges.size() + out_edges.size();
    }

    uint8_t in_degree() const {
        return in_edges.size();
    }

    uint8_t out_degree() const {
        return out_edges.size();
    }

    bool has_in_edge(id_t id) const {
        for (auto in_id : in_edges) {
            if (id == in_id) return true;
        }
        return false;
    }

    bool has_out_edge(id_t id) const {
        for (auto out_id : out_edges) {
            if (id == out_id) return true;
        }
        return false;
    }
};


typedef CompactNode * CompactNodePtr;
typedef DecisionNode * DecisionNodePtr;
typedef UnitigNode * UnitigNodePtr;

class cDBG : public KmerClient {

    std::unordered_map<hash_t, std::unique_ptr<DecisionNode>> decision_nodes;
    std::unordered_map<id_t, std::unique_ptr<UnitigNode>> unitig_nodes;
    std::unordered_map<hash_t, id_t> unitig_id_map;

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

    cDBG(uint16_t K) :
        KmerClient(K),
        _n_updates(0),
        _unitig_id_counter(UNITIG_START_ID),
        _n_unitig_nodes(0) {

    }

    DecisionNode* build_decision_node(hash_t hash, const string& kmer) {
        DecisionNode * dnode = get_decision_node(hash);
        if (dnode == nullptr) {
            unique_ptr<DecisionNode> dnode_ptr = make_unique<DecisionNode>(hash, kmer);
            decision_nodes.insert(make_pair(hash, std::move(dnode_ptr)));
            dnode = dnode_ptr.get();
        }
        return dnode;
    }

    DecisionNode* get_decision_node(hash_t hash) {
        auto search = decision_nodes.find(hash);
        if (search != decision_nodes.end()) {
            return search->second.get();
        }
        return nullptr;
    }

    template<class ShifterType>
    vector<DecisionNode*> get_decision_nodes(const string& sequence) {
        KmerIterator<ShifterType> iter(sequence, this->_K);
        vector<DecisionNode*> result;
        while(!iter.done()) {
            hash_t h = iter.next();
            DecisionNode * dnode;
            if ((dnode = get_decision_node(h)) != nullptr) {
                result.push_back(dnode);
            }
        }

        return result;
    }

    UnitigNode * build_unitig_node(HashVector& tags,
                                   const string& sequence,
                                   DecisionNode * left_dnode=nullptr,
                                   DecisionNode * right_dnode=nullptr) {

        unique_ptr<UnitigNode> unode = make_unique<UnitigNode>(_unitig_id_counter,
                                                               sequence);
        _unitig_id_counter++;
        _n_unitig_nodes++;

        unode->tags.insert(std::end(unode->tags), std::begin(tags), std::end(tags));
        unitig_nodes.insert(make_pair(unode->node_id, std::move(unode)));

        if (left_dnode != nullptr) {
            if (!left_dnode->has_out_edge(unode->node_id)) {
                left_dnode->out_edges.push_back(unode->node_id);
            }
        }
        if (right_dnode != nullptr) {
            if (!right_dnode->has_in_edge(unode->node_id)) {
                right_dnode->in_edges.push_back(unode->node_id);
            }
        }

        return unode.get();
    }
};


template <class GraphType>
class StreamingCompactor : public AssemblerMixin<GraphType> {

public:

    GraphType * dbg;
    cDBG cdbg;

    using AssemblerType = AssemblerMixin<GraphType>;
    using AssemblerType::seen;
    using AssemblerType::get_left;
    using AssemblerType::get_right;
    using AssemblerType::degree_left;
    using AssemblerType::degree_right;
    using AssemblerType::count_nodes;
    using AssemblerType::filter_nodes;

    StreamingCompactor(GraphType * dbg) :
        AssemblerMixin<GraphType>(dbg),
        dbg(dbg),
        cdbg(dbg->K()) {

    }

    string compactify(const string& seed) {
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

    bool is_decision_node(const string& node,
                          uint8_t& degree) {
        this->set_cursor(node);
        return is_decision_node(degree);
    }

    bool is_decision_node(const string& node) {
        this->set_cursor(node);
        return this->degree_left() > 1 || this->degree_right() > 1;
    }

    void find_decision_nodes(const string& sequence,
                             vector<uint32_t>& decision_positions,
                             HashVector& decision_hashes,
                             vector<NeighborBundle>& decision_neighbors) {

        KmerIterator<typename GraphType::shifter_type> iter(sequence, this->_K);
        size_t pos = 0;
        while(!iter.done()) {
            hash_t h = iter.next();
            NeighborBundle neighbors;
            if (get_decision_neighbors(iter.shifter,
                                       sequence.substr(pos, this->_K),
                                       neighbors)) {

                decision_neighbors.push_back(neighbors);
                decision_positions.push_back(pos);
                decision_hashes.push_back(h);
            }

            ++pos;
        }
    }

    bool get_decision_neighbors(typename GraphType::shifter_type& shifter,
                                const string& root_kmer,
                                NeighborBundle& result) {

        vector<shift_t> left_neighbors = filter_nodes(shifter.gather_left());
        vector<shift_t> right_neighbors = filter_nodes(shifter.gather_right());

        if (left_neighbors.size() > 1 || right_neighbors.size() > 1) {
            KmerVector left_kmers, right_kmers;
            for (auto neighbor : left_neighbors) {
                left_kmers.push_back(kmer_t(neighbor.hash,
                                            neighbor.symbol
                                            + root_kmer.substr(this->_K-1)));
            }
            for (auto neighbor : right_neighbors) {
                right_kmers.push_back(kmer_t(neighbor.hash,
                                             root_kmer.substr(1, this->_K-1)
                                             + neighbor.symbol));
                                      
            }
            result = std::make_pair(left_kmers, right_kmers);
            return true;
        } else {
            return false;
        }
    }

    void find_disturbed_decision_nodes(const string& sequence,
                                       vector<DecisionNode*>& disturbed_dnodes,
                                       vector<NeighborBundle>& disturbed_neighbors) {
        
        KmerIterator<typename GraphType::shifter_type> iter(sequence, this->_K);

        // first we have to search for induced decision nodes in-incident to
        // the first k-mer
        vector<shift_t> left_neighbors = filter_nodes(iter.shifter.gather_left());
        typename GraphType::shifter_type flank_shifter(iter.shifter);
        const char suffix = sequence[this->_K-1];
        for (shift_t left_neighbor : left_neighbors) {
            flank_shifter.shift_left(left_neighbor.symbol);
            NeighborBundle decision_neighbors;
            string kmer = flank_shifter.get_cursor();

            if (get_decision_neighbors(flank_shifter,
                                       kmer,
                                       decision_neighbors)) {
                DecisionNode * dnode;
                dnode = cdbg.build_decision_node(flank_shifter.get(), kmer);
                disturbed_dnodes.push_back(dnode);
                disturbed_neighbors.push_back(decision_neighbors);
            }
            flank_shifter.shift_right(suffix);
        }

        size_t pos = 0;
        while(!iter.done()) {
            hash_t h = iter.next();
            NeighborBundle decision_neighbors;
            string kmer = sequence.substr(pos, this->_K);
            if (get_decision_neighbors(iter.shifter,
                                       kmer,
                                       decision_neighbors)) {
                DecisionNode * dnode;
                dnode = cdbg.build_decision_node(h, kmer);
                disturbed_dnodes.push_back(dnode);
                disturbed_neighbors.push_back(decision_neighbors);
            }
            ++pos;
        }

        // and get the right flanking nodes
        vector<shift_t> right_neighbors = filter_nodes(iter.shifter.gather_right());
        const char prefix = sequence[sequence.length() - this->_K];
        for (shift_t right_neighbor : right_neighbors) {
            iter.shifter.shift_right(right_neighbor.symbol);
            NeighborBundle decision_neighbors;
            string kmer = iter.shifter.get_cursor();

            if (get_decision_neighbors(iter.shifter,
                                       kmer,
                                       decision_neighbors)) {
                DecisionNode * dnode;
                dnode = cdbg.build_decision_node(iter.shifter.get(), kmer);
                disturbed_dnodes.push_back(dnode);
                disturbed_neighbors.push_back(decision_neighbors);
            }
            iter.shifter.shift_left(prefix);
        }
    }

    bool insert_sequence(const string& sequence,
                         vector<uint32_t>& decision_positions,
                         HashVector& decision_hashes,
                         vector<NeighborBundle>& decision_neighbors) {

        if (!dbg->add_sequence(sequence)) {
            return false; // if there were no new k-mers, nothing to do
        } else {
            find_decision_nodes(sequence,
                                decision_positions,
                                decision_hashes,
                                decision_neighbors);
            return true;
        }
    }

    void update_linear(const string& sequence) {
    
    }

    /*
    void update(const string& sequence,
                vector<uint32_t>& decision_positions,
                HashVector& decision_hashes,
                vector<NeighborBundle>& decision_neighbors) {

        std::set<hash_t> updated_from;
        vector<kmer_t> update_left_queue;
        vector<kmer_t> update_right_queue;

        for (size_t node_idx=0; 
            node_idx<decision_positions.size(); ++node_idx) {
            hash_t node_hash = decision_hashes[node_idx];
            uint32_t kmer_pos = decision_positions[node_idx];
            const string kmer = sequence.substr(kmer_pos, this->_K);
            DecisionNode * decision_node = cdbg.build_decision_node(node_hash,
                                                                    kmer);
            NeighborBundle current_neighbors = decision_neighbors[node_idx];
            if (decision_node->left_degree() != 
                current_neighbors.first.size() ||
                decision_node->right_degree() != 
                current_neighbors.second.size()) {
                // Then this is either a new decision node or it 
                // got a new neighbor
                
                for (auto neighbor : current_neighbors.first) {
                    // left neighbors
                    if (new_kmers.count(neighbor.hash) && 
                        !started_from.count(neighbor.hash)) {
                        update_left(decision_node, neighbor, started_from);
                    }
                }

                for (auto neighbor : current_neighbors.second) {
                    // right neighbors
                    if (new_kmers.count(neighbor.hash) &&
                        !started_from.count(neighbor.hash)) {
                        update_left(decision_node, neighbor, started_from);
                    }
                }
            }
        }
    }
    */
};



}


#endif
