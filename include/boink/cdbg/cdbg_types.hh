/* cdbg/cdbg_types.hh -- cdbg node types
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_CDBG_TYPES_HH
#define BOINK_CDBG_TYPES_HH

#include <string>
#include <iostream>
#include <sstream>

#include "boink/boink.hh"
#include "oxli/kmer_hash.hh"

namespace boink {
namespace cdbg {

#define NULL_ID             ULLONG_MAX
#define UNITIG_START_ID     0

typedef uint64_t id_t;

enum node_meta_t {
    FULL,
    TIP,
    ISLAND,
    CIRCULAR,
    LOOP,
    TRIVIAL,
    DECISION
};


inline const char * node_meta_repr(node_meta_t meta) {
    switch(meta) {
        case FULL:
            return "FULL";
        case TIP:
            return "TIP";
        case ISLAND:
            return "ISLAND";
        case CIRCULAR:
            return "CIRCULAR";
        case TRIVIAL:
            return "TRIVIAL";
        case LOOP:
            return "LOOP";
        case DECISION:
            return "DECISION";
    }
}

enum cDBGFormat {
    GRAPHML,
    EDGELIST,
    FASTA,
    GFA1
};


inline const std::string cdbg_format_repr(cDBGFormat fmt) {
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


enum update_meta_t {
    BUILD_UNODE,
    BUILD_DNODE,
    DELETE_UNODE,
    SPLIT_UNODE,
    EXTEND_UNODE,
    CLIP_UNODE,
    MERGE_UNODES
};


class CompactNode {

protected:

    node_meta_t _meta;

public:

    const id_t node_id;
    id_t component_id;
    std::string sequence;
    
    CompactNode(id_t node_id,
                const std::string& sequence,
                node_meta_t meta)
        : _meta(meta),
          node_id(node_id),
          component_id(NULL_ID),
          sequence(sequence)
    {
    }

    std::string revcomp() const {
        return oxli::_revcomp(sequence);
    }

    size_t length() const {
        return sequence.length();
    }

    const node_meta_t meta() const {
        return _meta;
    }

    friend bool operator== (const CompactNode& lhs, const CompactNode& rhs) {
        return lhs.node_id == rhs.node_id;
    }

    std::string get_name() const {
        return std::string("NODE") + std::to_string(node_id);
    }

};


class DecisionNode: public CompactNode {

protected:

    bool _dirty;
    uint8_t _left_degree;
    uint8_t _right_degree;
    uint32_t _count;

public:

    DecisionNode(id_t node_id, const std::string& sequence)
        : CompactNode(node_id, sequence, DECISION),
          _dirty(true),
          _left_degree(0),
          _right_degree(0),
          _count(1)
    {    
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
      //<< " Dl=" << std::to_string(dn.left_degree())
      //<< " Dr=" << std::to_string(dn.right_degree())
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
               const std::string& sequence)
        : CompactNode(node_id, sequence, ISLAND),
          _left_end(left_end),
          _right_end(right_end) { 
    }

    void set_node_meta(node_meta_t new_meta) {
        _meta = new_meta;
    }

    const hash_t left_end() const {
        return _left_end;
    }

    void set_left_end(hash_t left_end) {
        _left_end = left_end;
    }

    void extend_right(hash_t right_end, const std::string& new_sequence) {
        sequence += new_sequence;
        _right_end = right_end;
    }

    void extend_left(hash_t left_end, const std::string& new_sequence) {
        sequence = new_sequence + sequence;
        _left_end = left_end;
    }

    const hash_t right_end() const {
        return _right_end;
    }

    void set_right_end(hash_t right_end) {
        _right_end = right_end;
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
      << " sequence=" << un.sequence
      << " length=" << un.sequence.length()
      << " meta=" << node_meta_repr(un.meta())
      << ">";
    return o;
}


typedef CompactNode * CompactNodePtr;
typedef DecisionNode * DecisionNodePtr;
typedef UnitigNode * UnitigNodePtr;



}
}
#endif
