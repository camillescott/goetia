/**
 * (c) Camille Scott, 2019
 * File   : cdbg.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 31.08.2019
 */


#ifndef GOETIA_CDBG_HH
#define GOETIA_CDBG_HH

#include <algorithm>
#include <cstdint>
#include <chrono>
#include <memory>
#include <mutex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>

// save diagnostic state
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wchar-subscripts"
#include "goetia/storage/sparsepp/spp.h"
#include "goetia/storage/phmap/phmap.h"
#pragma GCC diagnostic pop

#include "goetia/goetia.hh"
#include "goetia/metrics.hh"
#include "goetia/traversal/unitig_walker.hh"
#include "goetia/hashing/kmeriterator.hh"
#include "goetia/hashing/hashextender.hh"
#include "goetia/hashing/rollinghashshifter.hh"
#include "goetia/storage/storage_types.hh"
#include "goetia/cdbg/cdbg_types.hh"
#include "goetia/cdbg/metrics.hh"
#include "goetia/dbg.hh"


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


namespace goetia {

template <class T>
struct cDBG;


template <template <class, class> class GraphType,
          class StorageType,
          class ShifterType>
struct cDBG<GraphType<StorageType, ShifterType>> {

public:

    typedef GraphType<StorageType, ShifterType> graph_type;

    typedef typename graph_type::shifter_type   shifter_type;
    typedef HashExtender<shifter_type> extender_type;
    typedef typename shifter_type::alphabet     alphabet;
    typedef typename shifter_type::hash_type    hash_type;
  	typedef typename hash_type::value_type      value_type;
    typedef typename shifter_type::kmer_type    kmer_type;

  //typedef UnitigWalker<GraphType<StorageType, ShifterType>> walker_type;


    class CompactNode {

    protected:

        node_meta_t _meta;

    public:

        const id_t node_id;
        id_t component_id;
        std::string sequence;
        
        CompactNode(id_t node_id,
                    const std::string& sequence,
                    node_meta_t meta);

        std::string revcomp() const {
            return alphabet::reverse_complement(sequence);
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

        DecisionNode(id_t node_id, const std::string& sequence);

        static std::shared_ptr<DecisionNode> build(const DecisionNode& other) {
            return std::make_shared<DecisionNode>(other.node_id, other.sequence);
        }

        static std::shared_ptr<DecisionNode> build(const DecisionNode * other) {
            return std::make_shared<DecisionNode>(other->node_id, other->sequence);
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

        friend inline std::ostream& operator<<(std::ostream& o, const DecisionNode& dn) {

            o << "<DNode ID/hash=" << dn.node_id << " k-mer=" << dn.sequence
              //<< " Dl=" << std::to_string(dn.left_degree())
              //<< " Dr=" << std::to_string(dn.right_degree())
              << " count=" << dn.count()
              << " dirty=" << dn.is_dirty() << ">";
            return o;
        }
    };


    class UnitigNode : public CompactNode {

    protected:

        hash_type _left_end, _right_end;
        using CompactNode::_meta;

    public:

        using CompactNode::sequence;
        std::vector<hash_type> tags;

        UnitigNode(id_t node_id,
                   hash_type left_end,
                   hash_type right_end,
                   const std::string& sequence,
                   node_meta_t meta = ISLAND);

        static std::shared_ptr<UnitigNode> build(const UnitigNode& other) {
            return std::make_shared<UnitigNode>(other.node_id,
                                                other.left_end(),
                                                other.right_end(),
                                                other.sequence,
                                                other.meta());
        }

        static std::shared_ptr<UnitigNode> build(const UnitigNode * other) {
            return std::make_shared<UnitigNode>(other->node_id,
                                                other->left_end(),
                                                other->right_end(),
                                                other->sequence,
                                                other->meta());
        }

        void set_node_meta(node_meta_t new_meta) {
            _meta = new_meta;
        }

        const hash_type left_end() const {
            return _left_end;
        }

        void set_left_end(hash_type left_end) {
            _left_end = left_end;
        }

        void extend_right(hash_type right_end, const std::string& new_sequence) {
            sequence += new_sequence;
            _right_end = right_end;
        }

        void extend_left(hash_type left_end, const std::string& new_sequence) {
            sequence = new_sequence + sequence;
            _left_end = left_end;
        }

        const hash_type right_end() const {
            return _right_end;
        }

        void set_right_end(hash_type right_end) {
            _right_end = right_end;
        }

        std::string repr() const {
            std::ostringstream os;
            os << *this;
            return os.str();
        }

        friend inline std::ostream& operator<<(std::ostream& o, const UnitigNode& un) {
            o << "<UNode ID=" << un.node_id
              << " left_end=" << un.left_end()
              << " right_end=" << un.right_end()
              << " sequence=" << un.sequence
              << " length=" << un.sequence.length()
              << " meta=" << node_meta_repr(un.meta())
              << ">";
            return o;
        }

    };

    typedef CompactNode * CompactNodePtr;
    typedef DecisionNode * DecisionNodePtr;
    typedef UnitigNode * UnitigNodePtr;

    class Graph {

        /* Map of k-mer hash --> DecisionNode. DecisionNodes take
         * their k-mer hash value as their Node ID.
         */
        typedef phmap::parallel_flat_hash_map<value_type,
                                     std::unique_ptr<DecisionNode>> dnode_map_t;
        typedef typename dnode_map_t::const_iterator dnode_iter_t;

        /* Map of Node ID --> UnitigNode. This is a container
         * for the UnitigNodes' pointers; k-mer maps are stored elsewhere,
         * mapping k-mers to Node IDs.
         */
        typedef phmap::parallel_flat_hash_map<id_t,
                                     std::unique_ptr<UnitigNode>> unode_map_t;
        typedef typename unode_map_t::const_iterator unode_iter_t;

    protected:

        // The actual k-mer hash --> DNode map
        dnode_map_t decision_nodes;

        // The actual ID --> UNode map
        unode_map_t unitig_nodes;
        // The map from Unitig end k-mer hashes to UnitigNodes
        phmap::parallel_flat_hash_map<value_type, UnitigNode*> unitig_end_map;
        // The map from dBG k-mer tags to UnitigNodes
        phmap::parallel_flat_hash_map<value_type, UnitigNode*> unitig_tag_map;

        //std::mutex dnode_mutex;
        //std::mutex unode_mutex;
        std::mutex mutex;

        // Counts the number of cDBG updates so far
        uint64_t _n_updates;
        // Counter for generating UnitigNode IDs
        uint64_t _unitig_id_counter;
        // Current number of Unitigs
        uint64_t _n_unitig_nodes;

        id_t     component_id_counter;

    public:

        const uint16_t K;
        std::shared_ptr<graph_type> dbg;
        std::shared_ptr<cDBGMetrics> metrics;

        Graph(std::shared_ptr<graph_type> dbg,
              uint64_t minimizer_window_size=8);

        std::unique_lock<std::mutex> lock_nodes() {
            return std::unique_lock<std::mutex>(mutex);
        }

        /* Utility methods for iterating DNode and UNode
         * data structures. Note that these are not thread-safe
         * (the caller will need to lock).
         */

        typename dnode_map_t::const_iterator dnodes_begin() const {
            return decision_nodes.cbegin();
        }

        typename dnode_map_t::const_iterator dnodes_end() const {
            return decision_nodes.cend();
        }

        typename unode_map_t::const_iterator unodes_begin() const {
            return unitig_nodes.cbegin();
        }

        typename unode_map_t::const_iterator unodes_end() const {
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

        uint64_t n_unitig_ends() const {
            return unitig_end_map.size();
        }

        /* Node query methods: separate query mechanisms for
         * decision nodes and unitig nodes.
         */

        CompactNode* query_cnode(hash_type hash)  {
            CompactNode * node = query_unode_end(hash);
            if (node == nullptr) {
                node = query_dnode(hash);
            }
            return node;
        }

        DecisionNode* query_dnode(hash_type hash)  {

            auto search = decision_nodes.find(hash);
            if (search != decision_nodes.end()) {
                return search->second.get();
            }
            return nullptr;
        }

        std::vector<DecisionNode*> query_dnodes(const std::string& sequence)  {

            KmerIterator<shifter_type> kmers(sequence, this->K);
            std::vector<DecisionNode*> result;
            while(!kmers.done()) {
                hash_type h = kmers.next();
                DecisionNode * dnode;
                if ((dnode = query_dnode(h)) != nullptr) {
                    result.push_back(dnode);
                }
            }
            
            return result;
        }

        UnitigNode * query_unode_end(hash_type end_kmer)  {
            auto search = unitig_end_map.find(end_kmer);
            if (search != unitig_end_map.end()) {
                return search->second;
            }
            return nullptr;
        }

        UnitigNode * query_unode_tag(hash_type hash)  {
            auto search = unitig_tag_map.find(hash);
            if (search != unitig_tag_map.end()) {
                return search->second;
            }
            return nullptr;
        }

        UnitigNode * query_unode_id(id_t id)  {
            auto search = unitig_nodes.find(id);
            if (search != unitig_nodes.end()) {
                return search->second.get();
            }
            return nullptr;
        }

        bool has_dnode(hash_type hash)   {
            auto search = decision_nodes.find(hash);
            if (search != decision_nodes.end()) {
                return true;
            }
            return false;   
        }

        bool has_unode_end(hash_type end_kmer)  {
            return unitig_end_map.count(end_kmer) != 0;
        }

        CompactNode * find_rc_cnode(CompactNode * root)  {

            std::string  rc_seq  = alphabet::reverse_complement(root->sequence.substr(0, this->K));
            hash_type        rc_hash = dbg->hash(rc_seq);
            CompactNode * rc_node = query_cnode(rc_hash);

            return rc_node;
        }

        /* Neighbor-finding and traversal.
         *
         */

        auto find_dnode_neighbors(DecisionNode* dnode)
            -> std::pair<std::vector<CompactNode*>,
                         std::vector<CompactNode*>>;

        auto find_unode_neighbors(UnitigNode * unode)
            -> std::pair<DecisionNode*, DecisionNode*>;

        auto traverse_breadth_first(CompactNode* root)
            -> std::vector<CompactNode*>;

        auto find_connected_components()
            -> phmap::parallel_flat_hash_map<id_t, std::vector<id_t>>;

        /*
         * Graph Mutation
         */

        auto recompute_node_meta(UnitigNode * unode)
            -> node_meta_t;

        auto switch_unode_ends(hash_type old_unode_end,
                               hash_type new_unode_end)
            -> UnitigNode *;

        auto build_dnode(hash_type hash,
                         const std::string& kmer)
            -> DecisionNode*;

        auto build_unode(const std::string& sequence,
                          std::vector<hash_type>& tags,
                          hash_type left_end,
                          hash_type right_end)
            -> UnitigNode *;

        void clip_unode(bool      clip_from,
                        hash_type old_unode_end,
                        hash_type new_unode_end);

        void extend_unode(bool               ext_dir,
                          const std::string& new_sequence,
                          hash_type old_unode_end,
                          hash_type new_unode_end,
                          std::vector<hash_type>& new_tags);

        void split_unode(id_t node_id,
                         size_t split_at,
                         std::string split_kmer,
                         hash_type new_right_end,
                         hash_type new_left_end);

        void merge_unodes(const std::string& span_sequence,
                          size_t n_span_kmers,
                          hash_type left_end,
                          hash_type right_end,
                          std::vector<hash_type>& new_tags);

        void delete_unode(UnitigNode * unode) {

            if (unode != nullptr) {
                pdebug("Deleting " << *unode);
                id_t id = unode->node_id;
                metrics->decrement_cdbg_node(unode->meta());
                for (hash_type tag: unode->tags) {
                    unitig_tag_map.erase(tag);
                }
                unitig_end_map.erase(unode->left_end());
                unitig_end_map.erase(unode->right_end());

                unitig_nodes.erase(id);
                unode = nullptr;
                _n_unitig_nodes--;
                _n_updates++;
                metrics->n_unodes--;
                metrics->n_deletes++;
            }
        }
    
        void delete_unodes_from_tags(std::vector<hash_type>& tags) {

            for (auto tag: tags) {
                UnitigNode * unode = query_unode_tag(tag);
                if (unode != nullptr) {
                    delete_unode(unode);
                }
            }
        }

        void delete_dnode(DecisionNode * dnode){

            if (dnode != nullptr) {
                pdebug("Deleting " << *dnode);
                id_t id = dnode->node_id;
                metrics->n_dnodes--;
                metrics->n_deletes++;
                
                decision_nodes.erase(id);
                dnode = nullptr;

                _n_updates++;
            }
        }

        /*
         * File output
         */

        void validate(const std::string& filename) {

            std::ofstream out;
            out.open(filename);
            auto vdbg = dbg->clone();

            auto lock = lock_nodes();

            for (auto it = unitig_nodes.begin(); it != unitig_nodes.end(); ++it) {
                auto unode = it->second.get();
                auto counts = dbg->query_sequence(unode->sequence);
                vdbg->insert_sequence(unode->sequence);
                if (std::any_of(counts.begin(), counts.end(), 
                                [](count_t i){ return i == 0; })) {
                    out << unode->node_id << ";"
                        << unode->left_end() << ";"
                        << unode->right_end() << ";"
                        << unode->sequence << ";"
                        //<< counts
                        << std::endl;
                }
            }

            for (auto it = dnodes_begin(); it != dnodes_end(); ++it) {
                auto dnode = it->second.get();
                vdbg->insert_sequence(dnode->sequence);
            }

            out << "cdbg unique " << vdbg->n_unique() << std::endl;
            out << "reads unique " << dbg->n_unique() << std::endl;

            out.close();

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
                    throw GoetiaException("Invalid cDBG format.");
            };
        }

        void write_fasta(const std::string& filename)  {
            std::ofstream out;
            out.open(filename);
            write_fasta(out);
            out.close();
        }

        void write_fasta(std::ofstream& out)  {
            auto lock = lock_nodes();

            for (auto it = unitig_nodes.begin(); it != unitig_nodes.end(); ++it) {
                out << ">ID=" << it->first 
                    << " L=" << it->second->sequence.length()
                    << " type=" << node_meta_repr(it->second->meta())
                    << std::endl
                    << it->second->sequence
                    << std::endl;
            }
        }

        void write_gfa1(const std::string& filename)  {
            std::ofstream out;
            out.open(filename);
            write_gfa1(out);
            out.close();
        }

        void write_gfa1(std::ofstream& out);
   
        void write_graphml(const std::string& filename,
                           const std::string graph_name="cDBG") {

            std::ofstream out;
            out.open(filename);
            write_graphml(out, graph_name);
            out.close();
        }

        void write_graphml(std::ofstream& out,
                           const std::string graph_name="cDBG") {

            /*
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
            */

            auto lock = lock_nodes();

            //id_t edge_counter = 0;
            for (auto it = decision_nodes.begin(); it != decision_nodes.end(); ++it) {
                //out << "<node id=\"n" << it->first << "\"/>" << std::endl;
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

            //for (auto it = unitig_nodes.begin(); it != unitig_nodes.end(); ++it) {
            //    out << "<node id=\"" << it->first << "\"/>" << std::endl;
            //}

            //out << "</graph>" << std::endl;
            //out << "</graphml>" << std::endl;
        }
    };

    static auto compute_connected_component_metrics(std::shared_ptr<Graph> cdbg,
                                                    size_t                 sample_size = 10000)
        -> std::tuple<size_t, size_t, size_t, std::vector<size_t>>;

    static auto compute_unitig_fragmentation(std::shared_ptr<Graph> cdbg,
                                             std::vector<size_t>    bins)
        -> std::vector<size_t>;

};

extern template class goetia::cDBG<goetia::dBG<goetia::BitStorage, goetia::FwdLemireShifter>>;
extern template class goetia::cDBG<goetia::dBG<goetia::BitStorage, goetia::CanLemireShifter>>;

extern template class goetia::cDBG<goetia::dBG<goetia::PHMapStorage, goetia::FwdLemireShifter>>;
extern template class goetia::cDBG<goetia::dBG<goetia::PHMapStorage, goetia::CanLemireShifter>>;

extern template class goetia::cDBG<goetia::dBG<goetia::SparseppSetStorage, goetia::FwdLemireShifter>>;
extern template class goetia::cDBG<goetia::dBG<goetia::SparseppSetStorage, goetia::CanLemireShifter>>;

}


#undef pdebug
#endif
