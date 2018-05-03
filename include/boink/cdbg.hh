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
#include <memory>
#include <mutex>
#include <limits>
#include <iostream>
#include <sstream>
#include <unordered_map>

#include "oxli/alphabets.hh"
#include "gfakluge/src/gfakluge.hpp"

#include "boink/boink.hh"
#include "boink/hashing.hh"
#include "boink/minimizers.hh"

#include "boink/events.hh"
#include "boink/event_types.hh"

# ifdef DEBUG_CDBG
#   define pdebug(x) do { std::ostringstream stream; \
                          stream << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          std::cerr << stream.str(); \
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
using namespace boink::events;
using namespace boink::event_types;

#define NULL_ID             ULLONG_MAX
#define UNITIG_START_ID     0
#define NULL_JUNCTION       make_pair(0,0)

enum cDBGFormat {
    GRAPHML,
    EDGELIST,
    ADJMAT,
    FASTA,
    GFA1
};


inline const string cdbg_format_repr(cDBGFormat fmt) {
    switch(fmt) {
        case GRAPHML:
            return "graphml";
        case EDGELIST:
            return "edgelist";
        case ADJMAT:
            return "adjmat";
        case FASTA:
            return "fasta";
        case GFA1:
            return "gfa1";
        default:
            return "FORMAT";
    }
}


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

    string get_name() const {
        return string("NODE") + std::to_string(node_id);
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

class cDBG : public KmerClient,
             public EventListener {

public:

    typedef std::unordered_map<hash_t,
                               std::unique_ptr<DecisionNode>> dnode_map_t;
    typedef dnode_map_t::const_iterator dnode_iter_t;
    typedef std::unordered_map<id_t,
                               std::unique_ptr<UnitigNode>> unode_map_t;
    typedef unode_map_t::const_iterator unode_iter_t;

protected:

    dnode_map_t decision_nodes;
    std::mutex dnode_mutex;

    unode_map_t unitig_nodes;
    std::mutex unode_mutex;

    std::unordered_map<junction_t, id_t, junction_hash> unitig_junction_map;
    std::unordered_map<hash_t, id_t> unitig_tag_map;

    uint64_t _n_updates;
    uint64_t _unitig_id_counter;
    uint64_t _n_unitig_nodes;

public:

    node_meta_counter meta_counter;

    cDBG(uint16_t K)
        : KmerClient(K),
          EventListener("cDBG"),
          _n_updates(0),
          _unitig_id_counter(UNITIG_START_ID),
          _n_unitig_nodes(0)
    {
        this->msg_type_whitelist.insert(boink::event_types::MSG_ADD_DNODE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_ADD_UNODE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_DELETE_UNODE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_INCR_DNODE_COUNT);
    }

    std::unique_lock<std::mutex> lock_dnodes() {
        return std::unique_lock<std::mutex>(dnode_mutex);
    }

    std::unique_lock<std::mutex> lock_unodes() {
        return std::unique_lock<std::mutex>(unode_mutex);
    }

    virtual void handle_msg(shared_ptr<Event> event) {
        switch(event->msg_type) {
            case boink::event_types::MSG_ADD_DNODE:
                {
                    auto * data = static_cast<BuildDNodeEvent*>(event.get());
                    auto lock = lock_dnodes();
                    this->build_dnode(data->hash, data->kmer);
                }
                return;
            case boink::event_types::MSG_ADD_UNODE:
                {
                    auto * data = static_cast<BuildUNodeEvent*>(event.get());
                    auto lock = lock_unodes();
                    this->build_unode(data->tags, data->sequence,
                                      data->left, data->right);
                }
                return;
            case boink::event_types::MSG_DELETE_UNODE:
                {
                    auto * data = static_cast<DeleteUNodeEvent*>(event.get());
                    auto lock = lock_unodes();
                    this->delete_unode(get_unode_from_id(data->node_id));
                }
                return;
            case boink::event_types::MSG_INCR_DNODE_COUNT:
                return;
            default:
                return;
        }
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

        // quick and dirty
        delete_unode_from_tags(tags);

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

    void delete_unode_from_tags(HashVector& tags) {
        for (auto tag: tags) {
            UnitigNode * unode = get_unode(tag);
            if (unode != nullptr) {
                delete_unode(unode);
            }
        }
    }

    void write(const std::string& filename, cDBGFormat format) {
        std::ofstream out;
        out.open(filename);
        write(out, format);
        out.close();
    }

    void write(std::ofstream& out, cDBGFormat format) {
        switch (format) {
            case GRAPHML:
                write_graphml(out);
                break;
            case EDGELIST:
                write_edge_list(out);
                break;
            case ADJMAT:
                write_adj_matrix(out);
                break;
            case FASTA:
                write_fasta(out);
                break;
            case GFA1:
                write_gfa1(out);
                break;
            default:
                throw BoinkException("Invalid cDBG format.");
        };
    }

    void write_fasta(const string& filename) {
        std::ofstream out;
        out.open(filename);
        write_fasta(out);
        out.close();
    }

    void write_fasta(std::ofstream& out) {
        auto lock = lock_unodes();

        for (auto it = unitig_nodes.begin(); it != unitig_nodes.end(); ++it) {
            out << ">ID=" << it->first 
                << " L=" << it->second->sequence.length()
                << " type=" << node_meta_repr(it->second->meta)
                << std::endl
                << it->second->sequence
                << std::endl;
        }
    }

    void write_gfa1(const string& filename) {
        std::ofstream out;
        out.open(filename);
        write_gfa1(out);
        out.close();
    }

    void write_gfa1(std::ofstream& out) {
        auto lock2 = lock_dnodes();
        auto lock1 = lock_unodes();

        gfak::GFAKluge gfa;
        for (auto it = unitig_nodes.begin(); it != unitig_nodes.end(); ++it) {
            gfak::sequence_elem s;
            s.sequence = it->second->sequence;
            s.name = it->second->get_name();

            gfak::opt_elem ln_elem;
            ln_elem.key = "LN";
            ln_elem.val = std::to_string(it->second->sequence.length());
            ln_elem.type = "i";
            s.opt_fields.push_back(ln_elem);

            gfa.add_sequence(s);
        }

        for (auto it = decision_nodes.begin(); it != decision_nodes.end(); ++it) {
            gfak::sequence_elem s;
            s.sequence = it->second->sequence;
            s.name = it->second->get_name();

            gfak::opt_elem ln_elem;
            ln_elem.key = "LN";
            ln_elem.val = std::to_string(it->second->sequence.length());
            ln_elem.type = "i";
            s.opt_fields.push_back(ln_elem);

            gfa.add_sequence(s);
        }

        for (auto it = decision_nodes.begin(); it != decision_nodes.end(); ++it) {
            string root = it->second->get_name();
            for (auto junc : it->second->left_juncs) {
                UnitigNode * in_node = get_unode(junc);
                if (in_node == nullptr) continue;

                gfak::link_elem l;
                l.source_name = in_node->get_name();
                l.sink_name = root;
                l.source_orientation_forward = true;
                l.sink_orientation_forward = true;
                l.cigar = std::to_string(_K) + "M";

                string link_name = "LINK" + junction_to_string(junc);
                gfak::opt_elem id_elem;
                id_elem.key = "ID";
                id_elem.type = "Z";
                id_elem.val = link_name;
                l.opt_fields["ID"] = id_elem;

                gfa.add_link(in_node->get_name(), l);
            }
            for (auto junc : it->second->right_juncs) {
                UnitigNode * out_node = get_unode(junc);
                if (out_node == nullptr) continue;

                gfak::link_elem l;
                l.source_name = root;
                l.sink_name = out_node->get_name();
                l.source_orientation_forward = true;
                l.sink_orientation_forward = true;
                l.cigar = std::to_string(_K) + "M";

                string link_name = "LINK" + junction_to_string(junc);
                gfak::opt_elem id_elem;
                id_elem.key = "ID";
                id_elem.type = "Z";
                id_elem.val = link_name;
                l.opt_fields["ID"] = id_elem;
                gfa.add_link(root, l);
            }
        }

        out << gfa << std::endl;

    }

    void write_edge_list(const string& filename) {
        std::ofstream out;
        out.open(filename);
        write_edge_list(out);
        out.close();
    }

    void write_edge_list(std::ofstream& out) {
        auto lock1 = lock_unodes();
        auto lock2 = lock_dnodes();
    }

    void write_graphml(const string& filename,
                       string graph_name="cDBG") {
        std::ofstream out;
        out.open(filename);
        write_graphml(out, graph_name);
        out.close();
    }

    void write_graphml(std::ofstream& out,
                       const string graph_name="cDBG") {

        out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
               "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" "
               "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
               "xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns "
               "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">"
            << std::endl; // the header, open <graphml>
        out << "<graph id=\"" << graph_name 
            << "\" edgedefault=\"directed\" "
            << "parse.maxindegree=\"4\" parse.maxoutdegree=\"4\">"
            << std::endl; // open <graph>

        auto lock1 = lock_unodes();
        auto lock2 = lock_dnodes();

        id_t edge_counter = 0;
        for (auto it = decision_nodes.begin(); it != decision_nodes.end(); ++it) {
            out << "<node id=\"n" << it->first << "\"/>" << std::endl;
            for (auto junc : it->second->left_juncs) {
                id_t in_neighbor = unitig_junction_map[junc];
                out << "<edge id=\"e" << edge_counter 
                    << "\" source=\"n" << in_neighbor
                    << "\" target=\"n" << it->first << "\"/>"
                    << std::endl;
                edge_counter++;
            }
            for (auto junc : it->second->right_juncs) {
                id_t out_neighbor = unitig_junction_map[junc];
                out << "<edge id=\"e" << edge_counter 
                    << "\" source=\"n" << it->first
                    << "\" target=\"n" << out_neighbor << "\"/>"
                    << std::endl;
                edge_counter++;
            }
        }

        for (auto it = unitig_nodes.begin(); it != unitig_nodes.end(); ++it) {
            out << "<node id=\"" << it->first << "\"/>" << std::endl;
        }

        out << "</graph>" << std::endl;
        out << "</graphml>" << std::endl;
    }

    void write_adj_matrix(const string& filename) {
        std::ofstream out;
        out.open(filename);
        write_adj_matrix(out);
        out.close();
    }

    void write_adj_matrix(std::ofstream& out) {

        auto lock1 = lock_unodes();
        auto lock2 = lock_dnodes();

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
