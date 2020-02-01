/**
 * (c) Camille Scott, 2019
 * File   : cdbg.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 31.08.2019
 */


#ifndef BOINK_CDBG_HH
#define BOINK_CDBG_HH

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
#include "boink/storage/sparsepp/spp.h"
#include "boink/parsing/gfakluge/gfakluge.hpp"
#pragma GCC diagnostic pop
#include "boink/utils/stringutils.h"

#include "boink/boink.hh"
#include "boink/events.hh"
#include "boink/event_types.hh"
#include "boink/metrics.hh"
#include "boink/reporting/reporters.hh"
#include "boink/traversal.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/hashextender.hh"
#include "boink/hashing/rollinghashshifter.hh"
#include "boink/storage/storage_types.hh"
#include "boink/cdbg/cdbg_types.hh"
#include "boink/cdbg/metrics.hh"
#include "boink/dbg.hh"


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


namespace boink {
namespace cdbg {

template <class T>
struct cDBG;


template <template <class, class> class GraphType, class StorageType, class ShifterType>
struct cDBG<GraphType<StorageType, ShifterType>> {

public:

    typedef GraphType<StorageType, ShifterType> graph_type;

    typedef typename graph_type::shifter_type   shifter_type;
    typedef hashing::HashExtender<shifter_type> extender_type;
    typedef typename shifter_type::alphabet     alphabet;
    typedef typename shifter_type::hash_type    hash_type;
	typedef typename hash_type::value_type      value_type;
    typedef typename shifter_type::kmer_type    kmer_type;

    typedef dBGWalker<GraphType<StorageType, ShifterType>> walker_type;


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

        DecisionNode(id_t node_id, const std::string& sequence)
            : CompactNode(node_id, sequence, DECISION),
              _dirty(true),
              _left_degree(0),
              _right_degree(0),
              _count(1)
        {    
        }

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
                   node_meta_t meta = ISLAND)
            : CompactNode(node_id, sequence, meta),
              _left_end(left_end),
              _right_end(right_end) { 
        }

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

    class Graph : public events::EventNotifier {

        /* Map of k-mer hash --> DecisionNode. DecisionNodes take
         * their k-mer hash value as their Node ID.
         */
        typedef spp::sparse_hash_map<value_type,
                                     std::unique_ptr<DecisionNode>> dnode_map_t;
        typedef typename dnode_map_t::const_iterator dnode_iter_t;

        /* Map of Node ID --> UnitigNode. This is a container
         * for the UnitigNodes' pointers; k-mer maps are stored elsewhere,
         * mapping k-mers to Node IDs.
         */
        typedef spp::sparse_hash_map<id_t,
                                     std::unique_ptr<UnitigNode>> unode_map_t;
        typedef typename unode_map_t::const_iterator unode_iter_t;

    protected:

        // The actual k-mer hash --> DNode map
        dnode_map_t decision_nodes;

        // The actual ID --> UNode map
        unode_map_t unitig_nodes;
        // The map from Unitig end k-mer hashes to UnitigNodes
        spp::sparse_hash_map<value_type, UnitigNode*> unitig_end_map;
        // The map from dBG k-mer tags to UnitigNodes
        spp::sparse_hash_map<value_type, UnitigNode*> unitig_tag_map;

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
              uint64_t minimizer_window_size=8)
            : K(dbg->K),
              EventNotifier(),
              dbg(dbg),
              _n_updates(0),
              _unitig_id_counter(UNITIG_START_ID),
              _n_unitig_nodes(0),
              component_id_counter(0)
        {
            metrics = std::make_shared<cDBGMetrics>();
        }

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

            hashing::KmerIterator<shifter_type> kmers(sequence, this->K);
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

        std::pair<std::vector<CompactNode*>,
                  std::vector<CompactNode*>> find_dnode_neighbors(DecisionNode* dnode) {

            std::vector<CompactNode*> left;
            std::vector<CompactNode*> right;
            auto neighbors = dbg->neighbors(dnode->sequence);

            for (auto shift : neighbors.first) {
                CompactNode * node = query_cnode(shift.hash);
                if (node != nullptr) {
                    left.push_back(node);
                }
            }

            for (auto shift : neighbors.second) {
                CompactNode * node = query_cnode(shift.hash);
                if (node != nullptr) {
                    right.push_back(node);
                }
            }

            return std::make_pair(left, right);
        }

        std::pair<DecisionNode*, DecisionNode*> find_unode_neighbors(UnitigNode * unode) {

            DecisionNode * left = nullptr, * right = nullptr;

            dbg->set_cursor(unode->sequence.c_str());
            auto left_shifts = dbg->left_extensions();

            dbg->set_cursor(unode->sequence.c_str() + unode->sequence.size() - this->K);
            auto right_shifts = dbg->right_extensions();

            uint8_t n_left = 0;
            for (auto shift : left_shifts) {
                DecisionNode * dnode;
                if ((dnode = query_dnode(shift.hash)) != nullptr) {
                    left = dnode;
                    n_left++;
                    pdebug("Found left d-node: " << *dnode);
                }
            }

            uint8_t n_right = 0;
            for (auto shift : right_shifts) {
                DecisionNode * dnode;
                if ((dnode = query_dnode(shift.hash)) != nullptr) {
                    right = dnode;
                    n_right++;
                    pdebug("Found right d-node: " << *dnode);
                }
            }

            return std::make_pair(left, right);
        }

        std::vector<CompactNode*> traverse_breadth_first(CompactNode* root) {

            std::set<id_t> seen;
            std::vector<CompactNode*> node_q( {root} );
            std::vector<CompactNode*> result;

            while(!node_q.empty()) {
                root = node_q.back();
                node_q.pop_back();

                if (seen.count(root->node_id)) {
                    continue;
                } else {
                    result.push_back(root);
                }

                if (root->meta() == DECISION) {
                    auto neighbors = find_dnode_neighbors((DecisionNode*)root);
                    for (auto neighbor : neighbors.first) {
                        node_q.push_back(neighbor);
                    }
                    for (auto neighbor : neighbors.second) {
                        node_q.push_back(neighbor);
                    }
                } else {
                    auto neighbors = find_unode_neighbors((UnitigNode*)root);
                    if (neighbors.first != nullptr) {
                        node_q.push_back(neighbors.first);
                    }
                    if (neighbors.second != nullptr) {
                        node_q.push_back(neighbors.second);
                    }
                }

                seen.insert(root->node_id);
            }

            return result;
        }

        spp::sparse_hash_map<id_t, std::vector<id_t>> find_connected_components() {
            auto lock = this->lock_nodes();

            spp::sparse_hash_set<id_t>                    seen;
            spp::sparse_hash_map<id_t, std::vector<id_t>> components;

            for (auto unitig_it = unodes_begin(); unitig_it != unodes_end(); unitig_it++) {
                auto root = unitig_it->second.get();
                if (!seen.count(root->node_id)) {
                    auto              component_nodes = traverse_breadth_first(root);
                    std::vector<id_t> component_node_ids;
                    id_t              component_id = root->component_id;

                    if (component_id == NULL_ID) {
                        component_id = component_id_counter;
                        ++component_id_counter;
                        root->component_id = component_id;
                    }

                    for (auto node : component_nodes) {
                        node->component_id = component_id;
                        seen.insert(node->node_id);
                        component_node_ids.push_back(node->node_id);
                    }

                    components[component_id] = component_node_ids;
                }
            }

            return components;
        }

        /*
         * Graph Mutation
         */

        node_meta_t recompute_node_meta(UnitigNode * unode) {

            pdebug("Recompute node meta for " << unode->node_id);
            if (unode->sequence.size() == this->K) {
                return TRIVIAL;
            } else if (unode->left_end() == unode->right_end()) {
                return CIRCULAR;
            } else {
                auto neighbors = find_unode_neighbors(unode);
                if (neighbors.first == neighbors.second) {
                    if (neighbors.first == nullptr) {
                        return ISLAND;
                    } else {
                        return LOOP;
                    }
                } else {
                    if (neighbors.first == nullptr || neighbors.second == nullptr) {
                        return TIP;
                    } else {
                        return FULL;
                    }
                }
            }
        }

        UnitigNode * switch_unode_ends(hash_type old_unode_end,
                                       hash_type new_unode_end) {

            auto unode_end_it = unitig_end_map.find(old_unode_end);
            if (unode_end_it == unitig_end_map.end()) {
                return nullptr;
            }

            UnitigNode * unode = unode_end_it->second;
            unitig_end_map.erase(unode_end_it);
            unitig_end_map.insert(std::make_pair(new_unode_end, unode));

            pdebug("Swap " << old_unode_end << " to " << new_unode_end
                   << " for " << unode->node_id);

            return unode;
        }

        DecisionNode* build_dnode(hash_type hash,
                                  const std::string& kmer){
            /* Build a new DecisionNode; or, if the given k-mer hash
             * already has a DecisionNode, do nothing.
             */
            auto lock = lock_nodes();
            DecisionNode * dnode = query_dnode(hash);
            if (dnode == nullptr) {
                pdebug("BUILD_DNODE: " << hash << ", " << kmer);
                decision_nodes.emplace(hash,
                                       std::move(std::make_unique<DecisionNode>(hash, kmer)));
                // the memory location changes after the move; get a fresh address
                dnode = query_dnode(hash);
                notify_history_new(dnode->node_id,
                                   dnode->sequence,
                                   dnode->meta());
                metrics->n_dnodes++;
                pdebug("BUILD_DNODE complete: " << *dnode);
            } else {
                pdebug("BUILD_DNODE: d-node for " << hash << " already exists.");
                dnode->incr_count();
            }
            return dnode;
        }

        UnitigNode * build_unode(const std::string& sequence,
                                 std::vector<hash_type>& tags,
                                 hash_type left_end,
                                 hash_type right_end) {

            auto lock = lock_nodes();
            id_t id = _unitig_id_counter;
            
            // Transfer the UnitigNode's ownership to the map;
            // get its new memory address
            unitig_nodes.emplace(id, std::move(std::make_unique<UnitigNode>(id,
                                                                       left_end,
                                                                       right_end,
                                                                       sequence)));
            UnitigNode * unode_ptr = unitig_nodes[id].get();

            _unitig_id_counter++;
            _n_unitig_nodes++;
            _n_updates++;
            metrics->n_unodes++;

            // Link up its new tags
            unode_ptr->tags.insert(std::end(unode_ptr->tags),
                                   std::begin(tags),
                                   std::end(tags));
            for (auto tag: tags) {
                unitig_tag_map.insert(std::make_pair(tag, unode_ptr));
            }
            unitig_end_map.insert(std::make_pair(left_end, unode_ptr));
            unitig_end_map.insert(std::make_pair(right_end, unode_ptr));

            auto unode_meta = recompute_node_meta(unode_ptr);
            unode_ptr->set_node_meta(unode_meta);
            metrics->increment_cdbg_node(unode_meta);

            notify_history_new(id, unode_ptr->sequence, unode_ptr->meta());
            pdebug("BUILD_UNODE complete: " << *unode_ptr);

            return unode_ptr;
        }

        void clip_unode(bool      clip_from,
                        hash_type old_unode_end,
                        hash_type new_unode_end) {
    
            auto lock = lock_nodes();

            auto unode = switch_unode_ends(old_unode_end, new_unode_end);
            assert(unode != nullptr);
            pdebug("CLIP: " << *unode << " from " << (clip_from == hashing::DIR_LEFT ? std::string("LEFT") : std::string("RIGHT")) <<
                   " and swap " << old_unode_end << " to " << new_unode_end);

            if (unode->sequence.length() == this->K) {
                metrics->decrement_cdbg_node(unode->meta());
                delete_unode(unode);
                pdebug("CLIP complete: deleted null unode.");
            } else {
                metrics->n_clips++;
                if (clip_from == hashing::DIR_LEFT) {
                    unode->sequence = unode->sequence.substr(1);
                    unode->set_left_end(new_unode_end);

                    metrics->decrement_cdbg_node(unode->meta());
                    auto meta = recompute_node_meta(unode);
                    metrics->increment_cdbg_node(meta);
                    unode->set_node_meta(meta);

                    notify_history_clip(unode->node_id, unode->sequence, unode->meta());
                    pdebug("CLIP complete: " << *unode);
                } else {
                    unode->sequence = unode->sequence.substr(0, unode->sequence.length() - 1);
                    unode->set_right_end(new_unode_end);

                    metrics->decrement_cdbg_node(unode->meta());
                    auto meta = recompute_node_meta(unode);
                    metrics->increment_cdbg_node(meta);
                    unode->set_node_meta(meta);

                    notify_history_clip(unode->node_id, unode->sequence, unode->meta());
                    pdebug("CLIP complete: " << *unode);
                }
            }

            ++_n_updates;
        }

        void extend_unode(bool               ext_dir,
                          const std::string& new_sequence,
                          hash_type old_unode_end,
                          hash_type new_unode_end,
                          std::vector<hash_type>& new_tags) {

            auto lock = lock_nodes();

            auto unode = switch_unode_ends(old_unode_end, new_unode_end);
            if (unode->meta() == TRIVIAL) {
                unitig_end_map.insert(std::make_pair(old_unode_end, unode));
            }

            assert(unode != nullptr); 

            pdebug("EXTEND: from " << old_unode_end << " to " << new_unode_end
                   << (ext_dir == hashing::DIR_LEFT ? std::string(" to LEFT") : std::string(" to RIGHT"))
                   << " adding " << new_sequence << " to"
                   << std::endl << *unode);
            
            if (ext_dir == hashing::DIR_RIGHT) {
                unode->extend_right(new_unode_end, new_sequence);
            } else {
                unode->extend_left(new_unode_end, new_sequence);
            }

            std::copy(new_tags.begin(), new_tags.end(), std::back_inserter(unode->tags));
            for (auto tag: new_tags) {
                unitig_tag_map.insert(std::make_pair(tag, unode));
            }

            metrics->n_extends++;
            metrics->decrement_cdbg_node(unode->meta());
            auto meta = recompute_node_meta(unode);
            metrics->increment_cdbg_node(meta);
            unode->set_node_meta(meta);
            ++_n_updates;

            notify_history_extend(unode->node_id, unode->sequence, unode->meta());
            pdebug("EXTEND complete: " << *unode);
        }

        void split_unode(id_t node_id,
                         size_t split_at,
                         std::string split_kmer,
                         hash_type new_right_end,
                         hash_type new_left_end) {

            UnitigNode * unode;
            std::string right_unitig;
            hash_type right_unode_right_end;

            {
                auto lock = lock_nodes();

                unode = query_unode_id(node_id);
                assert(unode != nullptr);
                if (unode->meta() == CIRCULAR) {
                    pdebug("SPLIT: (CIRCULAR), flanking k-mers will become ends, " << 
                           new_left_end << " will be left_end, " << new_right_end <<
                           " will be right_end" << std::endl);

                    split_at = unode->sequence.find(split_kmer);
                    pdebug("Split k-mer found at " << split_at);
                    unode->sequence = unode->sequence.substr(split_at + 1) +
                                      unode->sequence.substr((this->K - 1), split_at);
                    switch_unode_ends(unode->left_end(), new_left_end);
                    unitig_end_map.insert(std::make_pair(new_right_end, unode));

                    unode->set_left_end(new_left_end);
                    unode->set_right_end(new_right_end);

                    metrics->n_splits++;
                    unode->set_node_meta(FULL);
                    metrics->decrement_cdbg_node(CIRCULAR);
                    metrics->increment_cdbg_node(FULL);
                    ++_n_updates;

                    notify_history_split_circular(unode->node_id, unode->sequence, unode->meta());
                    pdebug("SPLIT complete (CIRCULAR): " << *unode);
                    return;

                }
                pdebug("SPLIT: " << new_right_end << " left of root, "
                        << new_left_end << " right of root, at " << split_at
                        << std::endl << *unode);

                assert((split_at != 0) && (split_at != unode->sequence.size() - this->K));
                right_unitig = unode->sequence.substr(split_at + 1);

                // set the left unode right end to the new right end
                right_unode_right_end = unode->right_end();
                switch_unode_ends(unode->right_end(), new_right_end);
                unode->set_right_end(new_right_end);
                unode->sequence = unode->sequence.substr(0, split_at + this->K - 1);
                
                metrics->n_splits++;
                metrics->decrement_cdbg_node(unode->meta());
                auto meta = recompute_node_meta(unode);
                metrics->increment_cdbg_node(meta);
                unode->set_node_meta(meta);
                ++_n_updates;
            }

            std::vector<hash_type> tags; // TODO: WARNING: broken
            auto new_node = build_unode(right_unitig,
                                        tags,
                                        new_left_end,
                                        right_unode_right_end);

            notify_history_split(unode->node_id, unode->node_id, new_node->node_id,
                             unode->sequence, new_node->sequence,
                             unode->meta(), new_node->meta());
            pdebug("SPLIT complete: " << std::endl << *unode << std::endl << *new_node);

        }

        void merge_unodes(const std::string& span_sequence,
                          size_t n_span_kmers,
                          hash_type left_end,
                          hash_type right_end,
                          std::vector<hash_type>& new_tags) {
            /* span_sequence is the (K * 2) - 2 sequence connecting the two unitigs
             *
             */

            UnitigNode *left_unode, *right_unode;
            std::string right_sequence;
            hash_type new_right_end;

            {
                auto lock = lock_nodes();

                auto left_unode_it = unitig_end_map.find(left_end);
                if (left_unode_it == unitig_end_map.end()) {
                    return;
                }

                auto right_unode_it = unitig_end_map.find(right_end);
                if (right_unode_it == unitig_end_map.end()) {
                    return;
                }

                left_unode = left_unode_it->second;
                right_unode = right_unode_it->second;
            }

            auto rid = right_unode->node_id;

            if (left_unode->node_id == right_unode->node_id) {
                pdebug("MERGE: CIRCULAR! Creating circular unitig, span is " << span_sequence);
                metrics->decrement_cdbg_node(right_unode->meta());
                metrics->n_circular_merges++;
                std::string extend = span_sequence;
                //if (n_span_kmers < this->K - 1) {
                    pdebug("Overlap between merged sequence, trimming right " << n_span_kmers);
                    extend = span_sequence.substr(this->K-1, n_span_kmers);
                //} 
                extend_unode(hashing::DIR_RIGHT,
                             extend,
                             left_end, // this is left_unode's right_end
                             left_unode->left_end(),
                             new_tags);
            } else {

                pdebug("MERGE: " << left_end << " to " << right_end
                       << " with " << span_sequence 
                       << std::endl << *left_unode << std::endl << *right_unode);

                if (n_span_kmers < this->K - 1) {
                    size_t trim = (this->K - 1) - n_span_kmers;
                    pdebug("Overlap between merged sequence, trimming right " << trim);
                    right_sequence = right_unode->sequence.substr(trim);
                } else {
                    pdebug("No overlap, adding segment sequence, " << n_span_kmers);
                    right_sequence = span_sequence.substr(this->K - 1, n_span_kmers - this->K + 1)
                                                           + right_unode->sequence;
                }
                std::copy(right_unode->tags.begin(), right_unode->tags.end(),
                          std::back_inserter(new_tags));
                new_right_end = right_unode->right_end();

                delete_unode(right_unode);
                extend_unode(hashing::DIR_RIGHT,
                             right_sequence,
                             left_end,
                             new_right_end,
                             new_tags);
                metrics->n_merges++;

            }
            
            notify_history_merge(left_unode->node_id, rid,
                             left_unode->node_id,
                             left_unode->sequence,
                             left_unode->meta());
            pdebug("MERGE complete: " << *left_unode);
        }


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
         * Event notification
         */

        void notify_history_new(id_t id, std::string& sequence, node_meta_t meta) {

            auto event = std::make_shared<HistoryNewEvent>();
            event->id = id;
            event->sequence = sequence;
            event->meta = meta;
            this->notify(event);
        }

        void notify_history_merge(id_t lparent, id_t rparent, id_t child,
                                  std::string& sequence, node_meta_t meta) {

            auto event = std::make_shared<HistoryMergeEvent>();
            event->lparent = lparent;
            event->rparent = rparent;
            event->child = child;
            event->meta = meta;
            event->sequence = sequence;
            this->notify(event);
        }

        void notify_history_extend(id_t id, std::string& sequence, node_meta_t meta) {

            auto event = std::make_shared<HistoryExtendEvent>();
            event->id = id;
            event->sequence = sequence;
            event->meta = meta;
            this->notify(event);
        }

        void notify_history_clip(id_t id, std::string& sequence, node_meta_t meta) {
    
            auto event = std::make_shared<HistoryClipEvent>();
            event->id = id;
            event->sequence = sequence;
            event->meta = meta;
            this->notify(event);
        }

        void notify_history_split(id_t parent, id_t lchild, id_t rchild,
                                  std::string& lsequence, std::string& rsequence,
                                  node_meta_t lmeta, node_meta_t rmeta) {

            auto event = std::make_shared<HistorySplitEvent>();
            event->parent = parent;
            event->lchild = lchild;
            event->rchild = rchild;
            event->lsequence = lsequence;
            event->rsequence = rsequence;
            event->lmeta = lmeta;
            event->rmeta = rmeta;
            this->notify(event);
        }

        void notify_history_split_circular(id_t id, std::string& sequence, node_meta_t meta) {

            auto event = std::make_shared<HistorySplitCircularEvent>();
            event->id = id;
            event->sequence = sequence;
            event->meta = meta;
            this->notify(event);
        }


        /*
         * File output
         */

        void validate(const std::string& filename) {

            std::ofstream out;
            out.open(filename);

            auto lock = lock_nodes();

            for (auto it = unitig_nodes.begin(); it != unitig_nodes.end(); ++it) {
                auto unode = it->second.get();
                auto counts = dbg->query_sequence(unode->sequence);
                if (std::any_of(counts.begin(), counts.end(), 
                                [](storage::count_t i){ return i == 0; })) {
                    out << unode->node_id << ";"
                        << unode->left_end() << ";"
                        << unode->right_end() << ";"
                        << unode->sequence << ";"
                        //<< counts
                        << std::endl;
                }
            }

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
                    throw BoinkException("Invalid cDBG format.");
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

        void write_gfa1(std::ofstream& out) {
    
            auto lock = lock_nodes();

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

            auto get_link_name = [](id_t l, id_t r) { return std::string("LINK-") +
                                                             std::to_string(l) +
                                                             "-" +
                                                             std::to_string(r); };

            for (auto it = decision_nodes.begin(); it != decision_nodes.end(); ++it) {
                std::string root = it->second->get_name();
                auto neighbors = find_dnode_neighbors(it->second.get());

                for (auto in_node : neighbors.first) {

                    gfak::link_elem l;
                    l.source_name = in_node->get_name();
                    l.sink_name = root;
                    l.source_orientation_forward = true;
                    l.sink_orientation_forward = true;
                    l.cigar = std::to_string(K) + "M";

                    std::string link_name = get_link_name(in_node->node_id, it->second->node_id);
                    gfak::opt_elem id_elem;
                    id_elem.key = "ID";
                    id_elem.type = "Z";
                    id_elem.val = link_name;
                    l.opt_fields["ID"] = id_elem;

                    gfa.add_link(in_node->get_name(), l);
                }
                for (auto out_node : neighbors.second) {

                    gfak::link_elem l;
                    l.source_name = root;
                    l.sink_name = out_node->get_name();
                    l.source_orientation_forward = true;
                    l.sink_orientation_forward = true;
                    l.cigar = std::to_string(K) + "M";

                    std::string link_name = get_link_name(it->second->node_id, out_node->node_id);
                    gfak::opt_elem id_elem;
                    id_elem.key = "ID";
                    id_elem.type = "Z";
                    id_elem.val = link_name;
                    l.opt_fields["ID"] = id_elem;
                    gfa.add_link(root, l);
                }
            }
            out << gfa;
        }

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


    class ComponentReporter : public boink::reporting::SingleFileReporter {
    public:

        class Metrics {


            public:

                metrics::Gauge                         n_components;
                metrics::Gauge                         max_component_size;
                metrics::Gauge                         min_component_size;

                Metrics()
                    : 
                      n_components              {"size", "all_components"},
                      max_component_size        {"size", "max_component"},
                      min_component_size        {"size", "min_component"}
                {
                }
        };

    private:

        std::shared_ptr<Graph>            cdbg;

        uint64_t                          min_component;
        uint64_t                          max_component;

        // how large of a sample to take from the component size distribution
        size_t                            sample_size;
        metrics::ReservoirSample<size_t>  component_size_sample;

        std::unique_ptr<ComponentReporter::Metrics> metrics;

    public:

        ComponentReporter(std::shared_ptr<Graph>                 cdbg,
                          const std::string&                     filename,
                          size_t                                 sample_size = 10000)
            : SingleFileReporter       (filename, "cDBG::ComponentReporter"),
              cdbg                     (cdbg),
              min_component            (ULLONG_MAX),
              max_component            (0),
              sample_size              (sample_size),
              component_size_sample    (sample_size)
        {
            _cerr(this->THREAD_NAME << " reporting at MEDIUM interval.");
            this->msg_type_whitelist.insert(events::MSG_TIME_INTERVAL);
            _output_stream << "read_n,n_components,max_component,min_component,sample_size,component_size_sample" << std::endl;

            metrics = std::make_unique<ComponentReporter::Metrics>();
        }

        static std::shared_ptr<ComponentReporter> build(std::shared_ptr<Graph>                 cdbg,
                                                        const std::string&                     filename,
                                                        size_t                                 sample_size = 10000) {
            return std::make_shared<ComponentReporter>(cdbg,
                                                       filename,
                                                       sample_size);
        
        }

        virtual void handle_msg(std::shared_ptr<events::Event> event) {
             if (event->msg_type == events::MSG_TIME_INTERVAL) {
                auto _event = static_cast<events::TimeIntervalEvent*>(event.get());
                if (_event->level == events::TimeIntervalEvent::MEDIUM ||
                    _event->level == events::TimeIntervalEvent::END) {
                    
                    this->recompute_components();
                    _output_stream << _event->t << ","
                                   << component_size_sample.get_n_sampled() << ","
                                   << max_component << ","
                                   << min_component << ","
                                   << component_size_sample.get_sample_size() << ","
                                   << "\"" << repr(component_size_sample.get_result()) << "\""
                                   << std::endl;
                }
            }       
        }

        void recompute_components() {
            auto time_start = std::chrono::system_clock::now();

            component_size_sample.clear();
            auto components = cdbg->find_connected_components();
            for (auto id_comp_pair : components) {
                size_t component_size = id_comp_pair.second.size();
                component_size_sample.sample(component_size);
                max_component = (component_size > max_component) ? component_size : max_component;
                min_component = (component_size < min_component) ? component_size : min_component;
            }

            metrics->n_components.store(components.size());
            metrics->max_component_size.store(max_component);
            metrics->min_component_size.store(min_component);

            auto time_elapsed = std::chrono::system_clock::now() - time_start;
            _cerr("Finished recomputing components. Elapsed time: " <<
                  std::chrono::duration<double>(time_elapsed).count());


        }
    };

    class HistoryReporter : public reporting::SingleFileReporter {
        private:

        id_t _edge_id_counter;
        spp::sparse_hash_map<id_t, std::vector<std::string>> node_history;

        public:
        HistoryReporter(const std::string& filename)
            : SingleFileReporter(filename, "cDBG::HistoryReporter"),
              _edge_id_counter(0)
        {
            _cerr(this->THREAD_NAME << " reporting continuously.");

            this->msg_type_whitelist.insert(events::MSG_HISTORY_NEW);
            this->msg_type_whitelist.insert(events::MSG_HISTORY_SPLIT);
            this->msg_type_whitelist.insert(events::MSG_HISTORY_SPLIT_CIRCULAR);
            this->msg_type_whitelist.insert(events::MSG_HISTORY_MERGE);
            this->msg_type_whitelist.insert(events::MSG_HISTORY_EXTEND);
            this->msg_type_whitelist.insert(events::MSG_HISTORY_CLIP);
            this->msg_type_whitelist.insert(events::MSG_HISTORY_DELETE);

            _output_stream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
                              "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" "
                              "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
                              "xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns "
                              "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">"
                           << std::endl // the header, open <graphml>
                           << "<graph id=\"cDBG_History_DAG\" edgedefault=\"directed\">" << std::endl
                           << "<key id=\"op\" for=\"edge\" attr.name=\"op\" attr.type=\"string\"/>" << std::endl
                           << "<key id=\"seq\" for=\"node\" attr.name=\"seq\" attr.type=\"string\"/>" << std::endl
                           << "<key id=\"meta\" for=\"node\" attr.name=\"meta\" attr.type=\"string\"/>" << std::endl
                           << "<key id=\"node_id\" for=\"node\" attr.name=\"node_id\" attr.type=\"long\"/>"
                           << std::endl; // open <graph>
        }

        static std::shared_ptr<HistoryReporter> build(const std::string& filename) {
            return std::make_shared<HistoryReporter>(filename);
        }

        virtual void handle_exit() {
            _output_stream << "</graph>" << std::endl;
            _output_stream << "</graphml>" << std::endl;
        }

        void write_node(std::string id, id_t boink_id, std::string node_meta, std::string sequence) {
            _output_stream << "<node id=\"" << id << "\">" << std::endl
                           << "    <data key=\"seq\">" << sequence << "</data>" << std::endl
                           << "    <data key=\"meta\">" << node_meta << "</data>" << std::endl
                           << "    <data key=\"node_id\">" << boink_id << "</data>" << std::endl
                           << "</node>" << std::endl;
        }

        void write_edge(std::string src, std::string dst, std::string op) {
            auto id = _edge_id_counter++;
            _output_stream << "<edge id=\"" << id << "\" source=\"" 
                           << src << "\" target=\"" << dst << "\">" << std::endl
                           << "    <data key=\"op\">" << op << "</data>" << std::endl
                           << "</edge>" << std::endl;
        }

        std::string add_node_edit(id_t node_id, cdbg::node_meta_t meta, std::string sequence) {
            auto change_num = node_history[node_id].size();
            std::string id = std::to_string(node_id) + "_" + std::to_string(change_num);
            node_history[node_id].push_back(id);
            write_node(id, node_id, std::string(node_meta_repr(meta)), sequence);
            return id;
        }

        std::string add_new_node(id_t node_id, cdbg::node_meta_t meta, std::string sequence) {
            std::string id = std::to_string(node_id) + "_0";
            if (node_history.count(node_id) == 0) {
                node_history[node_id] = std::vector<std::string>{id};
                write_node(id, node_id, std::string(node_meta_repr(meta)), sequence);
            }
            return id;
        }

        virtual void handle_msg(std::shared_ptr<events::Event> event) {
            if (event->msg_type == events::MSG_HISTORY_NEW) {
                auto _event = static_cast<HistoryNewEvent*>(event.get());
                add_new_node(_event->id, _event->meta, _event->sequence);

            } else if (event->msg_type == events::MSG_HISTORY_SPLIT) {
                auto _event = static_cast<HistorySplitEvent*>(event.get());
                
                std::string parent_id = node_history[_event->parent].back();
                std::string lid, rid;
                if (_event->lchild == _event->parent) {
                    lid = add_node_edit(_event->lchild, _event->lmeta, _event->lsequence);
                    rid = add_new_node(_event->rchild, _event->rmeta, _event->rsequence);
                } else {
                    lid = add_new_node(_event->lchild, _event->lmeta, _event->lsequence);
                    rid = add_node_edit(_event->rchild, _event->rmeta, _event->rsequence);
                }
                write_edge(parent_id, lid, std::string("SPLIT"));
                write_edge(parent_id, rid, std::string("SPLIT"));

            } else if (event->msg_type == events::MSG_HISTORY_MERGE) {
                auto _event = static_cast<HistoryMergeEvent*>(event.get());

                std::string l_parent_id = node_history[_event->lparent].back();
                std::string r_parent_id = node_history[_event->rparent].back();
                std::string child_id = add_node_edit(_event->child, _event->meta, _event->sequence);
                
                write_edge(l_parent_id, child_id, std::string("MERGE"));
                write_edge(r_parent_id, child_id, std::string("MERGE"));

            } else if (event->msg_type == events::MSG_HISTORY_EXTEND) {
                auto _event = static_cast<HistoryExtendEvent*>(event.get());

                std::string src = node_history[_event->id].back();
                std::string dst = add_node_edit(_event->id, _event->meta, _event->sequence);
                write_edge(src, dst, std::string("EXTEND"));

            } else if (event->msg_type == events::MSG_HISTORY_CLIP) {
                auto _event = static_cast<HistoryClipEvent*>(event.get());

                std::string src = node_history[_event->id].back();
                std::string dst = add_node_edit(_event->id, _event->meta, _event->sequence);
                write_edge(src, dst, std::string("CLIP"));

            } else if (event->msg_type == events::MSG_HISTORY_SPLIT_CIRCULAR) {
                auto _event = static_cast<HistorySplitCircularEvent*>(event.get());

                std::string src = node_history[_event->id].back();
                std::string dst = add_node_edit(_event->id, _event->meta, _event->sequence);
                write_edge(src, dst, std::string("SPLIT_CIRCULAR"));
            }
        }
    };


    class UnitigReporter : public reporting::SingleFileReporter {
    private:

        std::shared_ptr<Graph>        cdbg;
        std::vector<size_t>           bins;

    public:

        UnitigReporter(std::shared_ptr<Graph> cdbg,
                       const std::string&     filename,
                       std::vector<size_t>    bins)
            : SingleFileReporter       (filename, "cDBG::UnitigReporter"),
              cdbg                     (cdbg),
              bins                     (bins)
        {
            _cerr(this->THREAD_NAME << " reporting at MEDIUM interval.");
            this->msg_type_whitelist.insert(events::MSG_TIME_INTERVAL);
            
            _output_stream << "read_n";

            for (size_t bin = 0; bin < bins.size() - 1; bin++) {
                _output_stream << ", " << bins[bin] << "-" << bins[bin+1];
            }
            _output_stream << ", " << bins.back() << "-Inf";

            _output_stream << std::endl;
        }

        static std::shared_ptr<UnitigReporter> build(std::shared_ptr<Graph> cdbg,
                                                     const std::string&     filename,
                                                     std::vector<size_t>    bins) {
            return std::make_shared<UnitigReporter>(cdbg, filename, bins);
        }

        virtual void handle_msg(std::shared_ptr<events::Event> event)  {
            if (event->msg_type == events::MSG_TIME_INTERVAL) {
                auto _event = static_cast<events::TimeIntervalEvent*>(event.get());
                if (_event->level == events::TimeIntervalEvent::MEDIUM ||
                    _event->level == events::TimeIntervalEvent::END) {
                    
                    auto bin_sums = this->compute_bins();
                    auto row      = utils::StringUtils::join(bin_sums, ", ");
                    _output_stream << _event->t << ","
                                   << row 
                                   << std::endl;
                }
            }       
        }

        std::vector<size_t> compute_bins() {
            auto time_start = std::chrono::system_clock::now();
            auto lock       = cdbg->lock_nodes();
            _cerr("Summing unitig length bins...");

            std::vector<size_t> bin_sums(bins.size(), 0);

            for (auto it = cdbg->unodes_begin(); it != cdbg->unodes_end(); ++it) {
                auto seq_len = it->second->sequence.length();
                for (size_t bin_num = 0; bin_num < bins.size() - 1; bin_num++) {
                    if (seq_len >= bins[bin_num] && seq_len < bins[bin_num+1]) {
                        bin_sums[bin_num] += seq_len;
                        break;
                    }
                }
                if (seq_len > bins.back()) {
                    bins[bins.size() - 1] += seq_len;
                }
            }

            auto time_elapsed = std::chrono::system_clock::now() - time_start;
            _cerr("Finished summing unitig length bins. Elapsed time: " <<
                  std::chrono::duration<double>(time_elapsed).count());

            return bin_sums;
        }
    };


    class Writer : public reporting::MultiFileReporter {
    protected:

        std::shared_ptr<Graph> cdbg;
        cdbg::cDBGFormat format;

    public:

        Writer(std::shared_ptr<Graph> cdbg,
                   cdbg::cDBGFormat format,
                   const std::string& output_prefix)
            : MultiFileReporter(output_prefix,
                                "cDBGWriter[" + cdbg_format_repr(format) + "]"),
              cdbg(cdbg),
              format(format)
        {
            _cerr(this->THREAD_NAME << " reporting at COARSE interval.");

            this->msg_type_whitelist.insert(events::MSG_TIME_INTERVAL);
        }

        static std::shared_ptr<Writer> build(std::shared_ptr<Graph> cdbg,
                                             cdbg::cDBGFormat format,
                                             const std::string& output_prefix) {
            return std::make_shared<Writer>(cdbg, format, output_prefix);
        }

        virtual void handle_msg(std::shared_ptr<events::Event> event) {
            if (event->msg_type == events::MSG_TIME_INTERVAL) {
                auto _event = static_cast<events::TimeIntervalEvent*>(event.get());
                if (_event->level == events::TimeIntervalEvent::COARSE ||
                    _event->level == events::TimeIntervalEvent::END) {

                    std::ofstream& stream = this->next_stream(_event->t,
                                                              cdbg_format_repr(format));
                    std::string&   filename = this->current_filename();

                    _cerr(this->THREAD_NAME << ", t=" << _event->t <<
                          ": write cDBG to " << filename);
                    cdbg->write(stream, format);
                }
            }
        }
    };


};




}

extern template class cdbg::cDBG<dBG<storage::BitStorage, hashing::FwdLemireShifter>>;
extern template class cdbg::cDBG<dBG<storage::BitStorage, hashing::CanLemireShifter>>;

extern template class cdbg::cDBG<dBG<storage::SparseppSetStorage, hashing::FwdLemireShifter>>;
extern template class cdbg::cDBG<dBG<storage::SparseppSetStorage, hashing::CanLemireShifter>>;

extern template class cdbg::cDBG<dBG<storage::ByteStorage, hashing::FwdLemireShifter>>;
extern template class cdbg::cDBG<dBG<storage::ByteStorage, hashing::CanLemireShifter>>;

extern template class cdbg::cDBG<dBG<storage::NibbleStorage, hashing::FwdLemireShifter>>;
extern template class cdbg::cDBG<dBG<storage::NibbleStorage, hashing::CanLemireShifter>>;

extern template class cdbg::cDBG<dBG<storage::QFStorage, hashing::FwdLemireShifter>>;
extern template class cdbg::cDBG<dBG<storage::QFStorage, hashing::CanLemireShifter>>;

}

#undef pdebug
#endif
