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
        pdebug("Build unode " << *(unode.get()));
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
};


template <class GraphType>
class StreamingCompactor : public AssemblerMixin<GraphType> {

protected:

    uint64_t _minimizer_window_size;

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

    using ShifterType = typename GraphType::shifter_type;

    StreamingCompactor(GraphType * dbg,
                       uint64_t minimizer_window_size=8) :
        AssemblerMixin<GraphType>(dbg),
        _minimizer_window_size(minimizer_window_size),
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

    void compactify_right(Path& path, InteriorMinimizer<hash_t>& minimizer) {
        this->seen.clear();
        this->seen.insert(this->get());
        minimizer.update(this->get());
        if (degree_left() > 1) return; // make sure we don't start at dnode
        
        shift_t next;
        while (get_right(next) // because get_right implicit checks rdegree
               && !this->seen.count(next.hash)) {
            path.push_back(next.symbol);
            this->seen.insert(next.hash);
            minimizer.update(next.hash);
            if (degree_left() > 1) break;
        }
    }

    void compactify_right(Path& path) {
        this->seen.clear();
        this->seen.insert(this->get());
        if (degree_left() > 1) return; // make sure we don't start at dnode
        
        shift_t next;
        while (get_right(next) // because get_right implicit checks rdegree
               && !this->seen.count(next.hash)) {
            path.push_back(next.symbol);
            this->seen.insert(next.hash);
            if (degree_left() > 1) break;
        }
    }

    void compactify_left(Path& path, InteriorMinimizer<hash_t> minimizer) {
        this->seen.clear();
        this->seen.insert(this->get());
        minimizer.update(this->get());
        if (degree_right() > 1) return; // don't start at dnode

        shift_t next;
        while (get_left(next) // because get_left implicitly checks ldegree
               && !this->seen.count(next.hash)) {
            path.push_front(next.symbol);
            this->seen.insert(next.hash);
            minimizer.update(next.hash);
            if (degree_right() > 1) break;
        }
    }

    void compactify_left(Path& path) {
        this->seen.clear();
        this->seen.insert(this->get());
        if (degree_right() > 1) return; // don't start at dnode

        shift_t next;
        while (get_left(next) // because get_left implicitly checks ldegree
               && !this->seen.count(next.hash)) {
            path.push_front(next.symbol);
            this->seen.insert(next.hash);
            if (degree_right() > 1) break;
        }
    }

    bool insert_sequence(const string& sequence,
                         vector<uint32_t>& decision_positions,
                         HashVector& decision_hashes,
                         vector<NeighborBundle>& decision_neighbors) {

        if (!dbg->add_sequence(sequence)) {
            return false; // if there were no new k-mers, nothing to do
        } else {
            find_decision_kmers(sequence,
                                decision_positions,
                                decision_hashes,
                                decision_neighbors);
            return true;
        }
    }

    void update(const string& sequence) {
        if(!dbg->add_sequence(sequence)) {
            return;
        } else {
            vector<DecisionNode*> disturbed_dnodes;
            vector<NeighborBundle> disturbed_neighbors;
            find_disturbed_dnodes(sequence,
                                  disturbed_dnodes,
                                  disturbed_neighbors);
            if (disturbed_dnodes.size() == 0) {
                _update_linear(sequence);
            } else {
                _update_from_dnodes(sequence,
                                    disturbed_dnodes,
                                    disturbed_neighbors);
            }
        }
    }

    void _update_linear(const string& sequence) {
        pdebug("Update linear");
        this->set_cursor(sequence.substr(0, this->_K));
        Path segment;
        segment.insert(segment.end(), sequence.begin(), sequence.end());
        this->compactify_left(segment);
        hash_t stopped_at = this->get();
        char after = *(segment.begin()+this->_K);
        junction_t left_junc = make_pair(stopped_at,
                                         this->shift_right(after));


        this->set_cursor(sequence.substr(sequence.length()-this->_K,
                                         this->_K));
        this->compactify_right(segment);
        stopped_at = this->get();
        after = *(segment.end()-(this->_K)-1);
        junction_t right_junc = make_pair(this->shift_left(after),
                               stopped_at);

        string segment_seq = this->to_string(segment);
        HashVector tags;
        WKMinimizer<ShifterType> minimizer(_minimizer_window_size,
                                           this->_K);
        tags = minimizer.get_minimizer_values(segment_seq);
        cdbg.build_unode(tags, segment_seq, left_junc, right_junc);
    }

    void _update_from_dnodes(const string& sequence,
                             vector<DecisionNode*>& disturbed_dnodes,
                             vector<NeighborBundle>& disturbed_neighbors) {

        std::set<hash_t> updated_dnodes;
        std::set<junction_t> updated_junctions;
        pdebug(disturbed_dnodes.size() << " on d-node queue, " <<
               disturbed_neighbors.size() << " on neighbor queue");

        while(!disturbed_dnodes.empty()) {
            DecisionNode * root_dnode = disturbed_dnodes.back();
            NeighborBundle root_neighbors = disturbed_neighbors.back();
            disturbed_dnodes.pop_back();
            disturbed_neighbors.pop_back();

            pdebug("updating from " << *root_dnode);

            if (updated_dnodes.count(root_dnode->node_id)) {
                continue;
            }

            for (kmer_t left_neighbor : root_neighbors.first) {
                junction_t right_junc = make_pair(left_neighbor.hash,
                                                  root_dnode->node_id);
				pdebug("left neighbor: " << left_neighbor.kmer 
                        << " junction: " << right_junc);
                if (updated_junctions.count(right_junc)) {
                    pdebug("Already updated from, continuing");
                    continue;
                }
                
                HashVector tags;
                std::string segment_seq;
                junction_t left_junc;
                if (cdbg.get_dnode(left_neighbor.hash) != nullptr) {
                    // trivial unode
                    left_junc = right_junc;
                    segment_seq = left_neighbor.kmer
                                  + root_dnode->sequence.back();
                    pdebug("Trivial u-node, junctions " << left_junc
                           << ", " << right_junc << ", sequence="
                           << segment_seq);
                } else {
                    this->set_cursor(left_neighbor.kmer);
                    Path this_segment;
                    this->get_cursor(this_segment);
                    pdebug("Start assembly left at " << this->get_cursor());
                    this_segment.push_back(root_dnode->sequence.back());

                    InteriorMinimizer<hash_t> minimizer(_minimizer_window_size);
                    this->compactify_left(this_segment, minimizer);
                    
                    hash_t stopped_at = this->get();
                    std::string after = std::string(this_segment.begin()+1,
                                                    this_segment.begin()+this->_K+1);
                    left_junc = make_pair(this->get(),
                                          this->shift_right(after.back()));
                    segment_seq = this->to_string(this_segment);
                    pdebug("Assembled left, stopped left of " << after
                            << ", shifting right on " << after.back()
                            << " with junction " << left_junc
                            << ", sequence=" << segment_seq);
                    tags = minimizer.get_minimizer_values();
                }
                updated_junctions.insert(left_junc);
                updated_junctions.insert(right_junc);

                cdbg.build_unode(tags, segment_seq, left_junc, right_junc);
            }

            for (kmer_t right_neighbor : root_neighbors.second) {
                junction_t left_junc = make_pair(root_dnode->node_id,
                                                 right_neighbor.hash);
				pdebug("right neighbor: " << right_neighbor.kmer 
                        << " junction: " << left_junc);
                if (updated_junctions.count(left_junc)) {
                    pdebug("Already updated from, continuing");
                    continue;
                }

                HashVector tags;
                std::string segment_seq;
                junction_t right_junc;
                if (cdbg.get_dnode(right_neighbor.hash) != nullptr) {
                    // trivial unode
                    right_junc = left_junc;
                    segment_seq = root_dnode->sequence.front()
                                  + right_neighbor.kmer;
                    pdebug("Trivial u-node, junctions " << left_junc
                           << ", " << right_junc << ", sequence="
                           << segment_seq);
                } else {
                    this->set_cursor(right_neighbor.kmer);
                    Path this_segment;
                    this->get_cursor(this_segment);
                    this_segment.push_front(root_dnode->sequence.front());

                    InteriorMinimizer<hash_t> minimizer(_minimizer_window_size);
                    this->compactify_right(this_segment, minimizer);

                    hash_t stopped_at = this->get();
                    std::string after = std::string(this_segment.end()-(this->_K)-1,
                                                    this_segment.end()-1);
                    right_junc = make_pair(this->shift_left(after.front()),
                                                      stopped_at);
                    segment_seq = this->to_string(this_segment);
                    pdebug("Assembled right, stopped right of " << after
                            << " with junction " << right_junc
                            << ", sequence=" << segment_seq);
                    tags = minimizer.get_minimizer_values();
                }
                updated_junctions.insert(left_junc);
                updated_junctions.insert(right_junc);

                cdbg.build_unode(tags, segment_seq, left_junc, right_junc);
            }

            updated_dnodes.insert(root_dnode->node_id);
        }
        pdebug("FINISHED _update_from_dnodes with " << updated_dnodes.size()
                << " updated from; " << cdbg.n_decision_nodes() 
                << " total d-nodes  and " << cdbg.n_unitig_nodes() << " u-nodes");
    }

    bool is_decision_kmer(uint8_t& degree) {
        // oops better put notes here
        uint8_t ldegree, rdegree;
        ldegree = this->degree_left();
        rdegree = this->degree_right();
        degree = ldegree + rdegree;
        return ldegree > 1 || rdegree > 1;
    }

    bool is_decision_kmer(const string& node,
                          uint8_t& degree) {
        this->set_cursor(node);
        return is_decision_kmer(degree);
    }

    bool is_decision_kmer(const string& node) {
        this->set_cursor(node);
        return this->degree_left() > 1 || this->degree_right() > 1;
    }

    void find_decision_kmers(const string& sequence,
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
                                            + root_kmer.substr(0, this->_K-1)));
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

    void find_disturbed_dnodes(const string& sequence,
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
                dnode = cdbg.build_dnode(flank_shifter.get(), kmer);
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
                pdebug("Found d-node " << h << ", " << kmer <<
                       " ldegree " << decision_neighbors.first.size() <<
                       " rdegree " << decision_neighbors.second.size());
                DecisionNode * dnode;
                dnode = cdbg.build_dnode(h, kmer);
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
                dnode = cdbg.build_dnode(iter.shifter.get(), kmer);
                disturbed_dnodes.push_back(dnode);
                disturbed_neighbors.push_back(decision_neighbors);
            }
            iter.shifter.shift_left(prefix);
        }
    }


};



}

#undef pdebug
#endif
