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
    FASTA,
    GFA1
};


inline const string cdbg_format_repr(cDBGFormat fmt) {
    switch(fmt) {
        case GRAPHML:
            return "graphml";
        case EDGELIST:
            return "edgelist";
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
    uint8_t _left_degree;
    uint8_t _right_degree;
    uint32_t _count;

public:

    DecisionNode(id_t node_id, const string& sequence) :
        CompactNode(node_id, sequence),
        _dirty(true),
        _left_degree(0),
        _right_degree(0),
        _count(1) {
        
    }

    const bool is_dirty() const {
        return _dirty;
    }

    void set_dirty(bool dirty) {
        _dirty = dirty;
    }

    const uint32_t count() const {
        return _count;
    }

    void incr_count() {
        _count++;
    }

    const uint8_t degree() const {
        return left_degree() + right_degree();
    }

    const uint8_t left_degree() const {
        return _left_degree;
    }

    void incr_left_degree() {
        _left_degree++;
    }

    const uint8_t right_degree() const {
        return _right_degree;
    }

    void incr_right_degree() {
        _right_degree++;
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
      << " count=" << dn.count()
      << " dirty=" << dn.is_dirty() << ">";
    return o;
}


class UnitigNode : public CompactNode {

protected:

    hash_t _left_end, _right_end;

public:

    HashVector tags;

    UnitigNode(id_t node_id,
               hash_t left_end,
               hash_t right_end,
               const string& sequence)
        : CompactNode(node_id, sequence),
          _left_end(left_end),
          _right_end(right_end) { 
    }

    const node_meta_t meta() const {
        // TODO: this should be deduced by the cDBG object
        return FULL;
    }

    const hash_t left_end() const {
        return _left_end;
    }

    void extend_right(hash_t right_end, const string& new_sequence) {
        sequence += new_sequence;
        _right_end = right_end;
    }

    void extend_left(hash_t left_end, const string& new_sequence) {
        sequence = new_sequence + sequence;
        _left_end = left_end;
    }

    const hash_t right_end() const {
        return _right_end;
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
      << " left_end=" << un.left_end()
      << " right_end=" << un.right_end()
      << " meta=" << node_meta_repr(un.meta())
      << ">";
    return o;
}


typedef CompactNode * CompactNodePtr;
typedef DecisionNode * DecisionNodePtr;
typedef UnitigNode * UnitigNodePtr;

class cDBG : public KmerClient {

public:

    /* Map of k-mer hash --> DecisionNode. DecisionNodes take
     * their k-mer hash value as their Node ID.
     */
    typedef std::unordered_map<hash_t,
                               std::unique_ptr<DecisionNode>> dnode_map_t;
    typedef dnode_map_t::const_iterator dnode_iter_t;

    /* Map of Node ID --> UnitigNode. This is a container
     * for the UnitigNodes' pointers; k-mer maps are stored elsewhere,
     * mapping k-mers to Node IDs.
     */
    typedef std::unordered_map<id_t,
                               std::unique_ptr<UnitigNode>> unode_map_t;
    typedef unode_map_t::const_iterator unode_iter_t;

protected:

    // The actual k-mer hash --> DNode map
    dnode_map_t decision_nodes;

    // The actual ID --> UNode map
    unode_map_t unitig_nodes;
    // The map from Unitig end k-mer hashes to UnitigNodes
    std::unordered_map<hash_t, UnitigNode*> unitig_end_map;
    // The map from dBG k-mer tags to UnitigNodes
    std::unordered_map<hash_t, UnitigNode*> unitig_tag_map;

    std::mutex dnode_mutex;
    std::mutex unode_mutex;

    // Counts the number of cDBG updates so far
    uint64_t _n_updates;
    // Counter for generating UnitigNode IDs
    uint64_t _unitig_id_counter;
    // Current number of Unitigs
    uint64_t _n_unitig_nodes;

public:

    // Container for cDBG metadata
    node_meta_counter meta_counter;

    cDBG(uint16_t K)
        : KmerClient(K),
          _n_updates(0),
          _unitig_id_counter(UNITIG_START_ID),
          _n_unitig_nodes(0) {}

    std::unique_lock<std::mutex> lock_dnodes() {
        return std::unique_lock<std::mutex>(dnode_mutex);
    }

    std::unique_lock<std::mutex> lock_unodes() {
        return std::unique_lock<std::mutex>(unode_mutex);
    }

    /* Utility methods for iterating DNode and UNode
     * data structures. Note that these are not thread-safe
     * (the caller will need to lock).
     */

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

    /* 
     * Accessor methods.
     */

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
        /* Build a new DecisionNode; or, if the given k-mer hash
         * already has a DecisionNode, do nothing.
         */
        auto lock = lock_dnodes();
        DecisionNode * dnode = query_dnode(hash);
        if (dnode == nullptr) {
            pdebug("Build d-node " << hash << ", " << kmer);
            unique_ptr<DecisionNode> dnode_ptr = make_unique<DecisionNode>(hash, kmer);
            decision_nodes.insert(make_pair(hash, std::move(dnode_ptr)));
            dnode = query_dnode(hash);
        }
        return dnode;
    }

    bool has_dnode(hash_t hash) {
        auto search = decision_nodes.find(hash);
        if (search != decision_nodes.end()) {
            return true;
        }
        return false;   
    }

    DecisionNode* query_dnode(hash_t hash) {
        auto search = decision_nodes.find(hash);
        if (search != decision_nodes.end()) {
            return search->second.get();
        }
        return nullptr;
    }

    template <class HashShifter>
    vector<DecisionNode*> query_dnodes(KmerIterator<HashShifter>& hash_iter) {
        vector<DecisionNode*> result;
        while(!hash_iter.done()) {
            hash_t h = hash_iter.next();
            DecisionNode * dnode;
            if ((dnode = query_dnode(h)) != nullptr) {
                result.push_back(dnode);
            }
        }
        
        return result;
    }

    UnitigNode * build_unode(const string& sequence,
                             HashVector& tags,
                             hash_t left_end,
                             hash_t right_end) {

        auto lock = lock_unodes();
        id_t id = _unitig_id_counter;
        unique_ptr<UnitigNode> unode = make_unique<UnitigNode>(id,
                                                               left_end,
                                                               right_end,
                                                               sequence);
        _unitig_id_counter++;
        _n_unitig_nodes++;

        // Transfer the UnitigNode's ownership to the map;
        // get its new memory address and return it
        unitig_nodes.insert(make_pair(id, std::move(unode)));
        UnitigNode * unode_ptr = unitig_nodes[id].get();

        // Link up its new tags
        unode_ptr->tags.insert(std::end(unode_ptr->tags),
                               std::begin(tags),
                               std::end(tags));
        for (auto tag: tags) {
            unitig_tag_map.insert(make_pair(tag, unode_ptr));
        }
        unitig_end_map.insert(make_pair(left_end, unode_ptr));
        unitig_end_map.insert(make_pair(right_end, unode_ptr));

        meta_counter.increment(unode_ptr->meta());

        pdebug("Built unode " << *unode_ptr);

        return unode_ptr;
    }

    void extend_unode(direction_t ext_dir,
                      const string& new_sequence,
                      hash_t old_unode_end,
                      hash_t new_unode_end,
                      HashVector& new_tags) {

        auto lock = lock_unodes();

        auto unode_end_it = unitig_end_map.find(old_unode_end);
        if (unode_end_it == unitig_end_map.end()) {
            return;
        }

        UnitigNode * unode = unode_end_it->second;
        unitig_end_map.erase(unode_end_it);
        unitig_end_map.insert(make_pair(new_unode_end, unode));

        if (ext_dir == DIR_RIGHT) {
            unode->extend_right(new_unode_end, new_sequence);
        } else {
            unode->extend_left(new_unode_end, new_sequence);
        }

        std::copy(new_tags.begin(), new_tags.end(), std::back_inserter(unode->tags));
        for (auto tag: new_tags) {
            unitig_tag_map.insert(make_pair(tag, unode));
        }

        ++_n_updates;
        pdebug("Extend unode " << *unode << " in dir " << ext_dir
               << " by " << new_sequence.size() << " bases, new end is "
               << new_unode_end << ", add " << new_tags.size() << " tags.");
    }

    void merge_unodes(const string& new_sequence,
                      hash_t left_end,
                      hash_t right_end,
                      HashVector& new_tags) {

        auto lock = lock_unodes();

        auto left_unode_it = unitig_end_map.find(left_end);
        if (left_unode_it == unitig_end_map.end()) {
            return;
        }

        auto right_unode_it = unitig_end_map.find(right_end);
        if (right_unode_it == unitig_end_map.end()) {
            return;
        }

        auto left_unode = left_unode_it->second;
        auto right_unode = right_unode_it->second;

        pdebug("Merge unodes " << *left_unode << " and " << *right_unode
               << " with delete and extend ops.");

        auto right_sequence = new_sequence + right_unode->sequence;
        std::copy(right_unode->tags.begin(), right_unode->tags.end(),
                  std::back_inserter(new_tags));
        hash_t new_right_end = right_unode->right_end();
        delete_unode(right_unode);

        extend_unode(DIR_RIGHT,
                     right_sequence,
                     left_end,
                     new_right_end,
                     new_tags);

        pdebug("Complete merge.");
    }

    UnitigNode * query_unode_end(hash_t end_kmer) {
        auto search = unitig_end_map.find(end_kmer);
        if (search != unitig_end_map.end()) {
            return search->second;
        }
        return nullptr;
    }

    bool has_unode_end(hash_t end_kmer) {
        return unitig_end_map.count(end_kmer) != 0;
    }

    UnitigNode * query_unode_tag(hash_t hash) {
        auto search = unitig_tag_map.find(hash);
        if (search != unitig_tag_map.end()) {
            return search->second;
        }
        return nullptr;
    }

    UnitigNode * query_unode_id(id_t id) {
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
            meta_counter.decrement(unode->meta());
            for (hash_t tag: unode->tags) {
                unitig_tag_map.erase(tag);
            }
            unitig_end_map.erase(unode->left_end());
            unitig_end_map.erase(unode->right_end());

            unitig_nodes.erase(id);
            unode = nullptr;
            _n_unitig_nodes--;
            _n_updates++;
        }
    }

    void delete_unodes_from_tags(HashVector& tags) {
        for (auto tag: tags) {
            UnitigNode * unode = query_unode_tag(tag);
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
                << " type=" << node_meta_repr(it->second->meta())
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

        /*
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
        */
        out << gfa << std::endl;

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
            /*
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
            */
        }

        for (auto it = unitig_nodes.begin(); it != unitig_nodes.end(); ++it) {
            out << "<node id=\"" << it->first << "\"/>" << std::endl;
        }

        out << "</graph>" << std::endl;
        out << "</graphml>" << std::endl;
    }



};

class AsyncCDBG : public cDBG,
                  public EventListener {

public:

    AsyncCDBG(uint16_t K)
        : cDBG(K),
          EventListener("cDBG")
    {
        /* Event types that the cDBG EventListener should filter for.
         * Other event types will be ignored by the EventListener
         * superclass on notify().
         */
        this->msg_type_whitelist.insert(boink::event_types::MSG_ADD_DNODE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_ADD_UNODE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_SPLIT_UNODE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_EXTEND_UNODE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_MERGE_UNODES);
        this->msg_type_whitelist.insert(boink::event_types::MSG_DELETE_UNODE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_INCR_DNODE_COUNT);
    }

    virtual void handle_msg(shared_ptr<Event> event) {
        /* Superclass override for handling the filtered messages.
         * For now, unrecognized messages are a no-op.
         * Follows the general pattern for the rest of boink:
         *  static_cast the Event to its subclassed type, lock
         *  underlying datastructures as appropriate, and handle
         *  updates.
         */
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
                    this->build_unode(data->sequence,
                                      data->tags,
                                      data->left_end,
                                      data->right_end);
                }
                return;
            case boink::event_types::MSG_DELETE_UNODE:
                {
                    auto * data = static_cast<DeleteUNodeEvent*>(event.get());
                    auto lock = lock_unodes();
                    this->delete_unode(query_unode_id(data->node_id));
                }
                return;
            case boink::event_types::MSG_INCR_DNODE_COUNT:
                return;
            default:
                return;
        }
    }

    void build_dnode_marker(hash_t hash) {
        // used to synchronously mark a d-node that
        // will be completely constructed eventually
        auto lock = lock_dnodes();
        if (query_dnode(hash) == nullptr) {
            decision_nodes.insert(make_pair(hash, nullptr));
        }
        pdebug("Built d-node marker: " << hash);
    }

    bool query_dnode_marker(hash_t hash) {
        // sync query for d-nodes that are either
        // constructed already or marked for construction
        // synchronously
        auto search = decision_nodes.find(hash);
        if (search != decision_nodes.end()) {
            return true;
        }
        return false;       
    }
};

}

#undef pdebug
#endif
