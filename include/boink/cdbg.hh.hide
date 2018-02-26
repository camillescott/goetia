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

#include "oxli/oxli.hh"
#include "oxli/kmer_hash.hh"
#include "oxli/hashtable.hh"
#include "oxli/hashgraph.hh"
#include "oxli/kmer_filters.hh"
#include "oxli/traversal.hh"
#include "oxli/assembler.hh"
#include "oxli/alphabets.hh"

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


class UnitigNodeFactory;
class UnitigNode {
    friend class UnitigNodeFactory;

public:

    string left_kmer;
    string right_kmer;
    uint32_t length;
    HashVector minimizers;
    node_meta_t meta;

    UnitigNode(string left_kmer, string right_kmer,
               uint32_t length) :
        left_kmer(left_kmer), right_kmer(right_kmer),
        length(length), meta(ISLAND), {}

    string rc_left() const {
        return _revcomp(left_kmer);
    }

    string rc_right() const {
        return _revcomp(right_kmer);
    }

    friend std::ostream& operator<<(std::ostream& stream,
                                     const UnitigNode& node) {
            stream << "<UnitigNode left=" << node.left_kmer
                   << " right=" << node.right_kmer
                   << " length=" << node.length
                   << " meta=" << node_meta_repr(node.meta)
                   << " n_minimizers=" << node.minimizers.size() << ">";
            return stream;
    }

};

typedef shared_ptr<UnitigNode> UnitigNodePtr;
typedef std::unordered_map<hash_t, UnitigNodePtr> HashUnitigMap;

typedef std::vector<UnitigNode> UnitigNodeVector;
typedef std::unordered_map<HashIntoType, UnitigNode*> TagEdgeMap;
typedef std::unordered_map<id_t, UnitigNode*> IDEdgeMap;
typedef std::pair<hash_t, UnitigNode*> TagEdgePair;
typedef std::set<TagEdgePair> TagEdgePairSet;
typedef std::set<UnitigNode*> UnitigNodeSet;

class DecisionNodeFactory;
class UnitigNodeFactory : public KmerFactory {
    friend class DecisionNodeFactory;
protected:

    uint64_t n_unitigs;
    uint64_t _n_updates;
    uint32_t _minimizer_w;

    HashUnitigMap unitigs;

public:

    UnitigNodeFactory(WordLength K) :

        KmerFactory(K), n_unitigs(0),
        _n_updates(0) {

        _minimizer_w = DEFAULT_TAG_DENSITY;    
    }

    uint64_t n_unitigs() const {
        return n_unitigs;
    }

    uint64_t n_updates() const {
        return _n_updates;
    }

    uint32_t minimizer_w() const {
        return _minimizer_w;
    }

    UnitigNode* build_edge(id_t left_id, id_t right_id,
                            compact_edge_meta_t edge_meta,
                            std::string edge_sequence) {

        UnitigNode* edge = new UnitigNode(left_id, right_id,
                                            _n_updates, edge_meta);
        compact_edges[_n_updates] = edge;

        pdebug("new compact edge: \n left=" << std::to_string(left_id) 
                << std::endl << " right=" << std::to_string(right_id)
                << std::endl << " meta=" << edge_meta_repr(edge_meta)
                << std::endl << " sequence   =" << edge_sequence
                << std::endl << " rc_sequence=" << _revcomp(edge_sequence)
                << std::endl << " start   =" << edge_sequence.substr(0, _ksize+1)
                << std::endl << " rc_start=" << _revcomp(edge_sequence.substr(0, _ksize+1))
                << std::endl << " end    =" 
                << edge_sequence.substr(edge_sequence.length()-_ksize-1, _ksize+1)
                << std::endl << " rc_end =" 
                << _revcomp(edge_sequence.substr(edge_sequence.length()-_ksize-1, _ksize+1)));

        edge->sequence = edge_sequence;
        n_compact_edges++;
        _n_updates++;
        return edge;
    }

    void link_tags_to_edge(UnitigNode * edge) {
        for (auto tag: edge->tags) {
            tags_to_edges[tag] = edge;
        }
    }

    UnitigNode* get_unitig_by_id(id_t id) {
        auto search = compact_edges.find(id);
        if (search != compact_edges.end()) {
            return search->second;
        }
        return nullptr;
    }

    void delete_edge(UnitigNode * edge) {
        //pdebug("attempt edge delete @" << edge);
        if (edge != nullptr) {
            pdebug("edge not null, proceeding");
            for (auto tag: edge->tags) {
                tags_to_edges.erase(tag);
            }
            compact_edges.erase(edge->edge_id);
            delete edge;
            edge = nullptr;
            n_compact_edges--;
            _n_updates++;
        }
    }

    void delete_edge_by_id(id_t id) {
        UnitigNode* e = get_unitig_by_id(id);
        delete_edge(e);
    }

    void delete_edge_by_tag(UHashSet& tags) {
        UnitigNode* edge = get_unitig(tags);
        delete_edge(edge);
    }

    void delete_edge_by_tag(HashIntoType tag) {
        UnitigNode* edge = get_unitig(tag);
        delete_edge(edge);
    }

    UnitigNode* get_unitig(HashIntoType tag) const {
        //pdebug("get compact edge from tag " << tag);
        auto search = tags_to_edges.find(tag);
        if (search != tags_to_edges.end()) {
            return search->second;
        }
        return nullptr;
    }

    bool get_tag_edge_pair(HashIntoType tag, TagEdgePair& pair) const {
        auto search = tags_to_edges.find(tag);
        if (search != tags_to_edges.end()) {
            pair = *search;
            return true;
        } else {
            return false;
        }
    }

    UnitigNode* get_unitig(UHashSet& tags) const {
        UnitigNode * edge = nullptr;
        for (auto tag: tags) {
            edge = get_unitig(tag);
            if (edge != nullptr) {
                break;
            }
        }
        return edge;
    }

    UnitigNode* get_unitig(KmerSet& kmers) const {
        UnitigNode * edge = nullptr;
        for (auto tag: kmers) {
            edge = get_unitig(tag);
            if (edge != nullptr) {
                break;
            }
        }
        return edge;
    }

    KmerFilter get_tag_stopper(TagEdgePair& te_pair,
                               bool& found_tag) { const
        KmerFilter stopper = [&] (const Kmer& node) {
            found_tag = get_tag_edge_pair(node, te_pair);
            return found_tag;
        };

        return stopper;
    }

    KmerHelper get_tag_gatherer(UHashSet& tags) { const

        KmerHelper gatherer = [&] (const Kmer& node) {
            UnitigNode * edge = get_unitig(node);
            if (edge != nullptr) {
                tags.insert(node);
            }
        };

        return gatherer;
    }

    void write_gml(const std::string filename,
                   const DecisionNodeFactory& nodes) const;
    void write_fasta(const std::string filename) const;

};


class DecisionNodeFactory;
class DecisionNode {
    friend class DecisionNodeFactory;
public:
    Kmer kmer;
    uint32_t count;
    const id_t node_id;
    std::string sequence;
    bool direction;

    UnitigNode* in_edges[4] = {nullptr, nullptr, nullptr, nullptr};
    UnitigNode* out_edges[4] = {nullptr, nullptr, nullptr, nullptr};

    DecisionNode(Kmer kmer, id_t node_id) : 
        kmer(kmer), count(0), node_id(node_id), direction(kmer.is_forward()) {}

    DecisionNode(Kmer kmer, std::string sequence, id_t node_id) : 
        kmer(kmer), count(0), sequence(sequence), node_id(node_id),
        direction(kmer.is_forward()) {}

    friend bool operator== (const DecisionNode& lhs, const DecisionNode& rhs) {
        return lhs.node_id == rhs.node_id;
    }

    std::string rc_sequence() const {
        return _revcomp(sequence);
    }

    bool delete_edge(UnitigNode* edge) {
        bool deleted = false;
        if (delete_in_edge(edge)) {
            deleted = true;
        }
        if (delete_out_edge(edge)) {
            deleted = true;
        }
        return deleted;
    }

    bool delete_in_edge(UnitigNode* edge) {
        for (uint8_t i=0; i<4; i++) {
            if (in_edges[i] == edge) {
                in_edges[i] = nullptr;
                return true;
            }
        }
        return false;
    }

    void add_in_edge(const char base, UnitigNode* edge) {
        //pdebug("add in edge to " << *this << ", base=" << base
        //        << ", edge: " << *edge);
        in_edges[twobit_repr(base)] = edge;
    }

    UnitigNode* get_in_edge(const char base) {
        return in_edges[twobit_repr(base)];
    }

    bool delete_out_edge(UnitigNode* edge) {
        for (uint8_t i=0; i<4; i++) {
            if (out_edges[i] == edge) {
                out_edges[i] = nullptr;
                return true;
            }
        }
        return false;
    }

    void add_out_edge(const char base, UnitigNode* edge) {
        //pdebug("add out edge to " << *this << ", base=" << base
        //        << ", edge: " << *edge);
        out_edges[twobit_repr(base)] = edge;
    }

    UnitigNode* get_out_edge(const char base) {
        return out_edges[twobit_repr(base)];
    }

    uint8_t degree() const {
        return out_degree() + in_degree();
    }

    uint8_t out_degree() const {
        uint8_t acc = 0;
        for (auto edge: out_edges) {
            if (edge != nullptr) {
                acc++;
            }
        }
        return acc;
    }

    uint8_t in_degree() const {
        uint8_t acc = 0;
        for (auto edge: in_edges) {
            if (edge != nullptr) {
                acc++;
            }
        }
        return acc;
    }

    friend std::ostream& operator<<(std::ostream& stream,
                                     const DecisionNode& node) {
            stream << "<DecisionNode ID=" << node.node_id << " Kmer=" << node.kmer.kmer_u
                   << " Sequence=" << node.sequence
                   << " rc_Sequence=" << node.rc_sequence()
                   << " Count=" << node.count << " in_degree=" 
                   << std::to_string(node.in_degree())
                   << " out_degree=" << std::to_string(node.out_degree()) << ">";
            return stream;
    }

    std::string edges_repr() {
        std::ostringstream os;
        os << *this << std::endl << "\tin_edges:" << std::endl;
        for (auto b : alphabets::DNA_SIMPLE) {
            UnitigNode* e = get_in_edge(b);
            if (e != nullptr) {
                os << "\t " << b << "=" << *e << std::endl;
            }
        }
        os << "\tout_edges:" << std::endl;
        for (auto b : alphabets::DNA_SIMPLE) {
            UnitigNode* e = get_out_edge(b);
            if (e != nullptr) {
                os << "\t " << b << "=" << *e << std::endl;
            }
        }
        return os.str();
    }
};

typedef std::vector<DecisionNode> DecisionNodeVector;

class DecisionNodeFactory : public KmerFactory {
    friend class UnitigNodeFactory;
protected:

    // map from HDN hashes to DecisionNode IDs
    HashIDMap kmer_id_map;
    // linear storage for DecisionNodes
    DecisionNodeVector compact_nodes;
    uint64_t n_compact_nodes;
    uint64_t _n_updates;

public:
    DecisionNodeFactory(WordLength K) : 
        KmerFactory(K), n_compact_nodes(0),
        _n_updates(0) {}

    uint64_t n_nodes() const {
        return n_compact_nodes;
    }
    
    uint64_t n_updates() const {
        return _n_updates;
    }

    // protected linear creation of DecisionNode
    // they should never be deleted, so this is straightforward
    DecisionNode* build_node(Kmer hdn) {
        pdebug("new compact node from " << hdn);
        DecisionNode * v = get_node_by_kmer(hdn);
        if (v == nullptr) {
            compact_nodes.emplace_back(hdn, n_compact_nodes);
            n_compact_nodes++;
            v = &(compact_nodes.back());
            v->sequence = _revhash(hdn, _ksize);
            kmer_id_map[hdn] = v->node_id;
            _n_updates++;
            pdebug("Allocate: " << *v);
        }
        return v;
    }

    DecisionNode* get_node_by_kmer(HashIntoType hdn) {
        auto search = kmer_id_map.find(hdn);
        if (search != kmer_id_map.end()) {
            id_t ID = search->second;
            return &(compact_nodes[ID]);
        }
        return nullptr;
    }

    DecisionNode* get_node_by_id(id_t id) {
        if (id >= compact_nodes.size()) {
            return nullptr;
        }
        return &(compact_nodes[id]);
    }

    DecisionNode* get_or_build_node(Kmer hdn) {
        DecisionNode* v = get_node_by_kmer(hdn);
        if (v != nullptr) {
            v->count += 1;
        } else {
            v = build_node(hdn);
            v->count = 1;
        }
        return v;
    }

    std::vector<DecisionNode*> get_nodes(const std::string& sequence) {
        //pdebug("get compact node IDs");
        KmerIterator kmers(sequence.c_str(), _ksize);
        std::vector<DecisionNode*> nodes;

        DecisionNode* node;

        while(!kmers.done()) {
            Kmer kmer = kmers.next();

            node = get_node_by_kmer(kmer);
            if (node != nullptr) {
                nodes.push_back(node);
            }
        }

        return nodes;
    }

    void unlink_edge(UnitigNode* edge) {
        pdebug("unlink edge " << *edge);
        DecisionNode *left, *right;
        left = get_node_by_id(edge->in_node_id);
        right = get_node_by_id(edge->out_node_id);
        if (left != nullptr) {
            // be lazy for now and use bidirectional delete
            left->delete_edge(edge);
            _n_updates++;
        }
        if (right != nullptr) {
            right->delete_edge(edge);
            _n_updates++;
        }
    }

    bool is_rc_from_left(DecisionNode* v, std::string& sequence) const {
        /* Check if sequence shares same canonical orientation with
         * v when coming from graph left, assuming sequence
         * does NOT include v.
         */
        const char * node_kmer = v->sequence.c_str();
        const char * _sequence = sequence.c_str();
        return strncmp(node_kmer, 
                       _sequence + sequence.size()-_ksize+1,
                       _ksize - 1) != 0;
    }

    bool get_pivot_from_left(DecisionNode* v,
                             std::string& sequence,
                             char& pivot_base) const {
        /* Check if sequence shared same canonical
         * orientation with v from graph left, assuming
         * sequence includes v
         */
        const char * node_kmer = v->sequence.c_str();
        const char * _segment = sequence.c_str();
        pivot_base = _segment[sequence.size()-_ksize-1];
        if (strncmp(node_kmer, 
                    _segment+sequence.size()-_ksize, 
                    _ksize-1) == 0) {
            // same canonical orientation
            return false;
        } else {
            // must have opposite canonical orientation
            pivot_base = complement(pivot_base);
            return true;
        }
    }

    bool add_edge_from_left(DecisionNode* v, UnitigNode* e) {
        char pivot_base;
        if (!get_pivot_from_left(v, e->sequence, pivot_base)) {
            // same canonical orientation
            pdebug("add in edge " << *e << " to node " << *v << " from " << pivot_base);
            v->add_in_edge(pivot_base, e);
            _n_updates++;
            return false;
        } else {
            // must have opposite canonical orientation
            pdebug("add out edge " << *e << " to node " << *v << " from " << pivot_base);
            v->add_out_edge(pivot_base, e);
            _n_updates++;
            return true;
        }
    }


    bool get_unitig_from_left(DecisionNode* v,
                            UnitigNode* &result_edge,
                            std::string& sequence) const {
        char pivot_base;
        if (!get_pivot_from_left(v, sequence, pivot_base)) {
            result_edge = v->get_in_edge(pivot_base);
            return false;
        } else {
            result_edge = v->get_out_edge(pivot_base);
            return true;
        }
    }

    bool is_rc_from_right(DecisionNode* v,
                          std::string& sequence) const {
        /* Check if sequence shared same canonical
         * orientation with v from graph right, assuming
         * sequence does NOT include v
         */
        const char * node_kmer = v->sequence.c_str();
        const char * _sequence = sequence.c_str();
        return strncmp(node_kmer+1, _sequence, _ksize-1) != 0;
    }

    bool get_pivot_from_right(DecisionNode* v,
                              std::string& sequence,
                              char& pivot_base) const {
        /* Find the "pivot base" between sequence and v
         * when sequence is from graph right, assuming
         * v contained in sequence
         */
        const char * node_kmer = v->sequence.c_str();
        const char * _segment = sequence.c_str();
        pivot_base = _segment[_ksize];
        if (strncmp(node_kmer+1, _segment+1, _ksize-1) == 0) {
            // same canonical orientation
            return false;
        } else {
            // must have opposite canonical orientation
            pivot_base = complement(pivot_base);
            return true;
        }
    }

    bool add_edge_from_right(DecisionNode* v, UnitigNode* e) {
        char pivot_base;
        if (!get_pivot_from_right(v, e->sequence, pivot_base)) {
            pdebug("add out edge " << *e << " to node " << *v << " from " << pivot_base);
            v->add_out_edge(pivot_base, e);
            _n_updates++;
            return false;
        } else {
            pdebug("add in edge " << *e << " to node " << *v << " from " << pivot_base);
            v->add_in_edge(pivot_base, e);
            _n_updates++;
            return true;
        }
    }

    bool get_unitig_from_right(DecisionNode* v,
                             UnitigNode* &result_edge,
                             std::string& sequence) const {
        char pivot_base;
        if (!get_pivot_from_right(v, sequence, pivot_base)) {
            result_edge = v->get_out_edge(pivot_base);
            return false;
        } else {
            result_edge = v->get_in_edge(pivot_base);
            return true;
        }

    }
};


class StreamingCompactor : public KmerFactory
{

protected:

    // map from tags to UnitigNodes
    DecisionNodeFactory nodes;
    UnitigNodeFactory edges;

    uint64_t n_sequences_added;

public:

    shared_ptr<Hashgraph> graph;
    
    StreamingCompactor(shared_ptr<Hashgraph> graph) :
        KmerFactory(graph->ksize()),
        nodes(graph->ksize()), edges(graph->ksize()),
        n_sequences_added(0), graph(graph)
    {
    }

    compact_edge_meta_t deduce_edge_meta(DecisionNode* in, DecisionNode* out) {
        compact_edge_meta_t edge_meta;
        if (in == nullptr && out == nullptr) {
           edge_meta = ISLAND;
        } else if ((out == nullptr) != (in == nullptr))  {
            edge_meta = TIP;
        } else {
            edge_meta = FULL;
        }
        return edge_meta;
    }

    uint64_t n_nodes() const {
        return nodes.n_nodes();
    }

    uint64_t n_edges() const {
        return edges.n_edges();
    }

    uint64_t n_updates() const {
        return nodes.n_updates() + edges.n_updates();
    }

    std::string report() const {
        std::ostringstream os;
        os << std::endl << "REPORT: StreamingCompactor(@" << this << " with "
           << "Hashgraph @" << graph.get() << ")" << std::endl;
        os << "  * " << n_nodes() << " cDBG nodes (HDNs)" << std::endl;
        os << "  * " << n_edges() << " cDBG edges" << std::endl;
        os << "  * " << n_sequences_added << " sequences added" << std::endl;
        return os.str();
    }


    DecisionNode* get_node_by_kmer(Kmer hdn) {
        return nodes.get_node_by_kmer(hdn);
    }

    DecisionNode* get_node_by_id(id_t id) {
        return nodes.get_node_by_id(id);
    }

    std::vector<DecisionNode*> get_nodes(const std::string& sequence) {
        return nodes.get_nodes(sequence);
    }

    UnitigNode* get_unitig(HashIntoType tag) const {
        return edges.get_unitig(tag);
    }

    bool get_tag_edge_pair(HashIntoType tag, TagEdgePair& pair) const {
        return edges.get_tag_edge_pair(tag, pair);
    }

    UnitigNode* get_unitig(UHashSet& tags) const {
        return edges.get_unitig(tags);
    }

    uint64_t consume_sequence(const std::string& sequence) {
        uint64_t prev_n_kmers = graph->n_unique_kmers();
        graph->consume_string(sequence);
        return graph->n_unique_kmers() - prev_n_kmers;
    }

    uint64_t consume_sequence_and_update(const std::string& sequence) {
        if (consume_sequence(sequence) > 0) {
            return update_compact_dbg(sequence);
        }
        return 0;
    }

    bool validate_segment(DecisionNode* root_node, DecisionNode* other_node,
                          UnitigNode* edge, std::string& sequence) {
        pdebug("validating " << *root_node << " with  " << *edge << ", " 
              << sequence << " and other node ID=" << 
              ((other_node != nullptr) ? other_node->node_id : NULL_ID));
        bool edge_valid = true;
        if (edge->meta == TIP) {
            if (other_node != nullptr) {
                edge_valid = false;
            }
            if (!((edge->in_node_id == root_node->node_id ||
                   edge->out_node_id == root_node->node_id) &&
                  edge->sequence.length() == sequence.length())) {
                edge_valid = false;
            }
        } else if (edge->meta == FULL) {
            if (other_node == nullptr) {
                edge_valid = false;
            } else {
                bool nodes_match;
                nodes_match = (edge->in_node_id == root_node->node_id && 
                               edge->out_node_id == other_node->node_id) ||
                              (edge->out_node_id == root_node->node_id &&
                               edge->in_node_id == other_node->node_id);
                if (!nodes_match) {
                    edge_valid = false;
                }
            }
        }
        pdebug("valid? = " << edge_valid);
        return edge_valid;
    }

    /* Update a compact dbg where there are no induced
     * HDNs
     */
    uint64_t update_compact_dbg_linear(const std::string& sequence,
                                       KmerSet& kmers) {
        pdebug("no induced HDNs, update linear...");
        uint64_t n_ops_before = n_updates();
        Kmer root_kmer = graph->build_kmer(sequence.substr(0, _ksize));

        UHashSet found_tags;
        KmerHelper tag_gatherer = edges.get_tag_gatherer(found_tags);
        CompactingAT<TRAVERSAL_LEFT> lcursor(graph.get(), root_kmer);
        lcursor.push_helper(tag_gatherer);
        CompactingAT<TRAVERSAL_RIGHT> rcursor(graph.get(), root_kmer);
        rcursor.push_helper(tag_gatherer);
        CompactingAssembler cassem(graph.get());

        std::string left_seq = cassem._assemble_directed(lcursor);
        std::string right_seq = cassem._assemble_directed(rcursor);
        std::string segment_seq = left_seq + right_seq.substr(_ksize);
        pdebug("linear segment: " << segment_seq);

        DecisionNode *left_node = nullptr, *right_node = nullptr;
        left_node = nodes.get_node_by_kmer(lcursor.cursor);
        right_node = nodes.get_node_by_kmer(rcursor.cursor);
        
        UnitigNode *left_edge = nullptr, *right_edge = nullptr;
        if (left_node != nullptr) {
            nodes.get_unitig_from_right(left_node, left_edge, segment_seq);
        }
        if (right_node != nullptr) {
            nodes.get_unitig_from_left(right_node, right_edge, segment_seq);
        }

        // meta for the induced edge
        compact_edge_meta_t edge_meta = deduce_edge_meta(left_node, 
                                                         right_node);

        if ((left_edge != nullptr && right_edge != nullptr &&
            left_edge->edge_id == right_edge->edge_id) 
            || left_edge != nullptr) {
            // account for a left edge OR a loop

            nodes.unlink_edge(left_edge);
            edges.delete_edge(left_edge);
        } else if (right_edge != nullptr) {
            nodes.unlink_edge(right_edge);
            edges.delete_edge(right_edge);
        }
        edges.delete_edge_by_tag(found_tags);

        if (edge_meta == ISLAND) {
            UnitigNode *new_edge = edges.build_edge(NULL_ID, NULL_ID,
                                                     edge_meta, segment_seq);
            // distribute tags
            uint32_t n_tags = ((segment_seq.length() - _ksize) 
                               / edges.tag_density()) + 1;
            uint32_t tag_dist = (segment_seq.length() - _ksize) / (n_tags + 1);
            uint32_t tag_pos = tag_dist;
            const char * _seq = segment_seq.c_str();
            pdebug(n_tags << " tags with dist " << tag_dist <<
                   " starting at " << tag_pos << ", seq length" << segment_seq.length());

            while (n_tags > 0) {
                new_edge->tags.insert(graph->hash_dna(_seq + tag_pos));
                tag_pos += tag_dist;
                n_tags--;
            }
            edges.link_tags_to_edge(new_edge);
            return n_updates() - n_ops_before;
        }

        id_t left_id, right_id;
        left_id = (left_node != nullptr) ? left_node->node_id : NULL_ID;
        right_id = (right_node != nullptr) ? right_node->node_id : NULL_ID;
        UnitigNode *new_edge = edges.build_edge(left_id, right_id,
                                                 edge_meta, segment_seq);
        if (left_node != nullptr) {
            nodes.add_edge_from_right(left_node, new_edge);
        }
        if (right_node != nullptr) {
            nodes.add_edge_from_left(right_node, new_edge);
        }
        
        return n_updates() - n_ops_before;
    }


    uint64_t update_compact_dbg(const std::string& sequence) {
        pdebug("update cDBG from " << sequence);
        n_sequences_added++;
        uint64_t n_ops_before = n_updates();

        // first gather up all k-mers that could have been disturbed --
        // k-mers in the read, and the neighbors of the flanking nodes
        KmerIterator kmers(sequence.c_str(), _ksize);
        KmerQueue disturbed_kmers;
        Kmer kmer = kmers.next();
        CompactingAT<TRAVERSAL_LEFT> lcursor(graph.get(), kmer);
        lcursor.neighbors(disturbed_kmers);
        while(!kmers.done()) {
            kmer = kmers.next();
            disturbed_kmers.push_back(kmer);
        }
        CompactingAT<TRAVERSAL_RIGHT> rcursor(graph.get(), kmer);
        rcursor.neighbors(disturbed_kmers);
        
        pdebug(disturbed_kmers.size() << " k-mers disturbed" << std::endl);
        
                // find the induced HDNs in the disturbed k-mers
        KmerSet induced_hdns;
        KmerSet disturbed_hdns;
        for (Kmer kmer : disturbed_kmers) {
            uint8_t l_degree, r_degree;
            l_degree = lcursor.degree(kmer);
            r_degree = rcursor.degree(kmer);
            if(l_degree > 1 || r_degree > 1) {
                pdebug("found HDN... " << kmer);
                DecisionNode* hdn = nodes.get_or_build_node(kmer);
                if (hdn->count == 1) { // just created
                    induced_hdns.insert(kmer);
                } else if (hdn->degree() != (l_degree + r_degree)) {
                    induced_hdns.insert(kmer);
                } else {
                    disturbed_hdns.insert(kmer);
                }
            }
        }
        pdebug(induced_hdns.size() << " induced HDNs");

        /* If there are no induced HDNs, we must have extended
         * a tip or merged two tips into a linear segment */
        if (induced_hdns.size() == 0 && disturbed_hdns.size() == 0) {
            return update_compact_dbg_linear(sequence, disturbed_kmers);
        } else if (induced_hdns.size() == 0) {
            induced_hdns.insert(disturbed_hdns.begin(), disturbed_hdns.end());
        }

        /* Update from all induced HDNs
         */
        CompactingAssembler cassem(graph.get());
        KmerQueue neighbors;
        while(!induced_hdns.empty()) {
            Kmer root_kmer = *induced_hdns.begin();
            induced_hdns.erase(root_kmer);

            DecisionNode* root_node = nodes.get_node_by_kmer(root_kmer);
            char root_front = root_node->sequence.front();
            char root_back = root_node->sequence.back();
            pdebug("searching from induced HDN: " << root_node->edges_repr());

            // check left (in) edges
            lcursor.neighbors(root_kmer, neighbors);
            pdebug("checking " << neighbors.size() << " left neighbors");
            while(!neighbors.empty()) {
                Kmer neighbor = neighbors.back();
                neighbors.pop_back();
                lcursor.cursor = neighbor;

                TagEdgePair tag_pair;
                bool found_tag = false;

                lcursor.push_filter(edges.get_tag_stopper(tag_pair, found_tag));
                std::string segment_seq = cassem._assemble_directed(lcursor);
                if (nodes.is_rc_from_left(root_node, segment_seq)) {
                    segment_seq = segment_seq  + complement(root_front);
                } else {
                    segment_seq = segment_seq + root_back;
                }
                pdebug("assembled segment: " << segment_seq << " length: " << 
                       segment_seq.length());

                // first check for a segment going this direction from root
                UnitigNode* segment_edge = nullptr;
                nodes.get_unitig_from_left(root_node, segment_edge, segment_seq);

                DecisionNode* left_node = nodes.get_node_by_kmer(lcursor.cursor);
                UnitigNode* left_out_edge = nullptr;
                if (left_node != nullptr) {
                    pdebug("found existing left node: " << *left_node);
                    nodes.get_unitig_from_right(left_node, left_out_edge, segment_seq);
                }

                // validate edge leaving root if it exists
                if (segment_edge != nullptr && left_out_edge != nullptr) {
                    pdebug("found edges leaving root and left node");
                    
                    if (segment_edge == left_out_edge && 
                        validate_segment(root_node, left_node, 
                                         segment_edge, segment_seq)) {
                        continue;
                    } else {
                        nodes.unlink_edge(segment_edge);
                        nodes.unlink_edge(left_out_edge);
                        edges.delete_edge(segment_edge);
                        edges.delete_edge(left_out_edge);
                    }
                } else if (left_out_edge != nullptr) {
                    // there was no edge from root, must be bad
                    pdebug("edge from left invalid, delete");
                    nodes.unlink_edge(left_out_edge);
                    edges.delete_edge(left_out_edge);
                } else if (segment_edge != nullptr) {
                    pdebug("found end leaving root node");
                    if (validate_segment(root_node, left_node,
                                         segment_edge, segment_seq)) {
                        continue;
                    } else {
                        pdebug("edge from root invalid, delete");
                        nodes.unlink_edge(segment_edge);
                        edges.delete_edge(segment_edge);
                    }
                }
                
                /*
                 * Should also keep a set of pair<Kmer,Kmer> to track resolved
                 * segments
                 */

                // not needed until tags used again
                //segment_seq = cassem._assemble_directed(lcursor) +
                //              segment_seq.substr(_ksize);

                // construct the compact edge
                compact_edge_meta_t edge_meta = (left_node == nullptr) 
                                                ? TIP : FULL;
                edge_meta = (segment_seq.length() == _ksize + 1 && edge_meta == FULL)
                            ? TRIVIAL : edge_meta;

                if (edge_meta == FULL || edge_meta == TRIVIAL) {
                    segment_edge = edges.build_edge(left_node->node_id, 
                                                    root_node->node_id,
                                                    edge_meta, 
                                                    segment_seq);
                    nodes.add_edge_from_right(left_node, segment_edge);
                } else {
                    segment_edge = edges.build_edge(NULL_ID, 
                                                    root_node->node_id,
                                                    edge_meta, 
                                                    segment_seq);
                }

                nodes.add_edge_from_left(root_node, segment_edge);
            }

            // now the right neighbors...
            rcursor.neighbors(root_kmer, neighbors);
            pdebug("checking " << neighbors.size() << " right neighbors");
            while(!neighbors.empty()) {
                Kmer neighbor = neighbors.back();
                neighbors.pop_back();
                rcursor.cursor = neighbor;
                pdebug("right neighbor: " << neighbor.repr(_ksize));

                TagEdgePair tag_pair;
                bool found_tag = false;

                rcursor.push_filter(edges.get_tag_stopper(tag_pair, found_tag));
                std::string segment_seq = cassem._assemble_directed(rcursor);
                if (nodes.is_rc_from_right(root_node, segment_seq)) {
                    segment_seq = complement(root_back) + segment_seq;
                } else {
                    segment_seq = root_front + segment_seq;
                }
                pdebug("assembled segment: " << segment_seq << " length: " << 
                       segment_seq.length());
                // first check for a segment going this direction from root
                UnitigNode* segment_edge = nullptr;
                nodes.get_unitig_from_right(root_node, segment_edge, segment_seq);

                DecisionNode* right_node = nodes.get_node_by_kmer(rcursor.cursor);
                UnitigNode* right_in_edge = nullptr;
                if (right_node != nullptr) {
                    nodes.get_unitig_from_left(right_node, right_in_edge, segment_seq);
                }

                // validate edge leaving root if it exists
                if (segment_edge != nullptr && right_in_edge != nullptr) {

                    
                    if (segment_edge == right_in_edge && 
                        validate_segment(root_node, right_node, 
                                         segment_edge, segment_seq)) {
                        continue;
                    } else {
                        nodes.unlink_edge(segment_edge);
                        nodes.unlink_edge(right_in_edge);
                        edges.delete_edge(segment_edge);
                        edges.delete_edge(right_in_edge);
                    }
                } else if (right_in_edge != nullptr) {
                    // there was no edge from root, must be bad
                    pdebug("edge from left invalid, delete");
                    nodes.unlink_edge(right_in_edge);
                    edges.delete_edge(right_in_edge);
                } else if (segment_edge != nullptr) {
                    if (validate_segment(root_node, right_node,
                                         segment_edge, segment_seq)) {
                        continue;
                    } else {
                        pdebug("edge from root invalid, delete");
                        nodes.unlink_edge(segment_edge);
                        edges.delete_edge(segment_edge);
                    }
                }

                compact_edge_meta_t edge_meta = (right_node == nullptr) ?
                                                  TIP : FULL;
                edge_meta = (segment_seq.length() == _ksize + 1 && edge_meta == FULL)
                            ? TRIVIAL : edge_meta;

                if (edge_meta == FULL || edge_meta == TRIVIAL) {
                    segment_edge = edges.build_edge(root_node->node_id, 
                                                    right_node->node_id,
                                                    edge_meta, 
                                                    segment_seq);
                    nodes.add_edge_from_left(right_node, segment_edge);
                } else {
                    segment_edge = edges.build_edge(root_node->node_id,
                                                    NULL_ID,
                                                    edge_meta, 
                                                    segment_seq);
                }

                nodes.add_edge_from_right(root_node, segment_edge);
            }

        }

        return n_updates() - n_ops_before;

    } // update_compact_dbg

    void write_gml(const std::string filename) const {
        edges.write_gml(filename, nodes);
    }

    void write_fasta(const std::string filename) const {
        edges.write_fasta(filename);
    }

};



}


#endif
