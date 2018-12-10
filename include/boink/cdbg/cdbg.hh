/* cdbg.hh -- compact de Bruijn Graph
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_CDBG_HH
#define BOINK_CDBG_HH

#include <algorithm>
#include <cstdint>
#include <memory>
#include <mutex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <tuple>

#include "gfakluge/gfakluge.hpp"
#include "sparsepp/spp.h"

#include "boink/assembly.hh"
#include "boink/boink.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/minimizers.hh"
#include "boink/storage/storage.hh"

#include "boink/events.hh"
#include "boink/event_types.hh"

#include "boink/cdbg/cdbg_types.hh"
#include "boink/cdbg/metrics.hh"

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

#define ifmetrics(x) do { if (metrics != nullptr) { x; } } while (0)

namespace boink {
namespace cdbg {

using boink::hashing::hash_t;
using boink::hashing::shift_t;
using boink::hashing::kmer_t;
using boink::hashing::KmerIterator;

using std::unique_ptr;
using std::make_unique;
using std::vector;
using std::pair;


template <class GraphType>
class cDBG : public hashing::KmerClient,
             public events::EventNotifier {

protected:

    using ShifterType = typename GraphType::shifter_type;
    using CompactorType = CompactorMixin<GraphType>;
    using MinimizerType = WKMinimizer<ShifterType>;

public:

    typedef GraphType graph_type;
    typedef ShifterType shifter_type;
    typedef CompactorType compactor_type;
    typedef MinimizerType minimizer_type;

    /* Map of k-mer hash --> DecisionNode. DecisionNodes take
     * their k-mer hash value as their Node ID.
     */
    typedef spp::sparse_hash_map<hash_t,
                                 std::unique_ptr<DecisionNode>> dnode_map_t;
    typedef dnode_map_t::const_iterator dnode_iter_t;

    /* Map of Node ID --> UnitigNode. This is a container
     * for the UnitigNodes' pointers; k-mer maps are stored elsewhere,
     * mapping k-mers to Node IDs.
     */
    typedef spp::sparse_hash_map<id_t,
                                 std::unique_ptr<UnitigNode>> unode_map_t;
    typedef unode_map_t::const_iterator unode_iter_t;

protected:

    shared_ptr<GraphType> dbg;

    // The actual k-mer hash --> DNode map
    dnode_map_t decision_nodes;

    // The actual ID --> UNode map
    unode_map_t unitig_nodes;
    // The map from Unitig end k-mer hashes to UnitigNodes
    spp::sparse_hash_map<hash_t, UnitigNode*> unitig_end_map;
    // The map from dBG k-mer tags to UnitigNodes
    spp::sparse_hash_map<hash_t, UnitigNode*> unitig_tag_map;

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

    shared_ptr<prometheus::Registry> pr_registry;

public:

    shared_ptr<cDBGMetrics> metrics;

    cDBG(shared_ptr<GraphType> dbg,
         shared_ptr<prometheus::Registry> metrics_registry=nullptr,
         uint64_t minimizer_window_size=8)
        : KmerClient(dbg->K()),
          EventNotifier(),
          dbg(dbg),
          _n_updates(0),
          _unitig_id_counter(UNITIG_START_ID),
          _n_unitig_nodes(0),
          component_id_counter(0),
          pr_registry(metrics_registry)
    {
        if(metrics_registry != nullptr) {
            pr_registry = metrics_registry;
            metrics = make_shared<cDBGMetrics>(pr_registry);
        }
    }

    std::unique_lock<std::mutex> lock_nodes() {
        return std::unique_lock<std::mutex>(mutex);
    }

    /* Utility methods for iterating DNode and UNode
     * data structures. Note that these are not thread-safe
     * (the caller will need to lock).
     */

    GraphType * get_dbg() {
        return dbg;
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

    CompactNode* query_cnode(hash_t hash) {
        CompactNode * node = query_unode_end(hash);
        if (node == nullptr) {
            node = query_dnode(hash);
        }
        return node;
    }

    DecisionNode* query_dnode(hash_t hash) {

        auto search = decision_nodes.find(hash);
        if (search != decision_nodes.end()) {
            return search->second.get();
        }
        return nullptr;
    }

    vector<DecisionNode*> query_dnodes(const std::string& sequence) {

        KmerIterator<ShifterType> kmers(sequence, this->_K);
        vector<DecisionNode*> result;
        while(!kmers.done()) {
            hash_t h = kmers.next();
            DecisionNode * dnode;
            if ((dnode = query_dnode(h)) != nullptr) {
                result.push_back(dnode);
            }
        }
        
        return result;
    }

    UnitigNode * query_unode_end(hash_t end_kmer) {
        auto search = unitig_end_map.find(end_kmer);
        if (search != unitig_end_map.end()) {
            return search->second;
        }
        return nullptr;
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

    bool has_dnode(hash_t hash) {
        auto search = decision_nodes.find(hash);
        if (search != decision_nodes.end()) {
            return true;
        }
        return false;   
    }

    bool has_unode_end(hash_t end_kmer) {
        return unitig_end_map.count(end_kmer) != 0;
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

        return make_pair(left, right);
    }

    std::pair<DecisionNode*, DecisionNode*> find_unode_neighbors(UnitigNode * unode) {
        DecisionNode * left = nullptr, * right = nullptr;
        CompactorType compactor(dbg);

        compactor.set_cursor(unode->sequence.c_str());
        auto left_shifts = compactor.gather_left();

        compactor.set_cursor(unode->sequence.c_str() + unode->sequence.size() - this->_K);
        auto right_shifts = compactor.gather_right();

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

    vector<CompactNode*> traverse_breadth_first(CompactNode* root) {
        std::set<id_t> seen;
        vector<CompactNode*> node_q( {root} );
        vector<CompactNode*> result;

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

    node_meta_t recompute_node_meta(UnitigNode * unode) {
        pdebug("Recompute node meta for " << unode->node_id);
        if (unode->sequence.size() == this->_K) {
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

    UnitigNode * switch_unode_ends(hash_t old_unode_end,
                                   hash_t new_unode_end) {

        auto unode_end_it = unitig_end_map.find(old_unode_end);
        if (unode_end_it == unitig_end_map.end()) {
            return nullptr;
        }

        UnitigNode * unode = unode_end_it->second;
        unitig_end_map.erase(unode_end_it);
        unitig_end_map.insert(make_pair(new_unode_end, unode));

        pdebug("Swap " << old_unode_end << " to " << new_unode_end
               << " for " << unode->node_id);

        return unode;
    }

    DecisionNode* build_dnode(hash_t hash,
                              const std::string& kmer) {
        /* Build a new DecisionNode; or, if the given k-mer hash
         * already has a DecisionNode, do nothing.
         */
        auto lock = lock_nodes();
        DecisionNode * dnode = query_dnode(hash);
        if (dnode == nullptr) {
            pdebug("BUILD_DNODE: " << hash << ", " << kmer);
            decision_nodes.emplace(hash,
                                   std::move(make_unique<DecisionNode>(hash, kmer)));
            // the memory location changes after the move; get a fresh address
            dnode = query_dnode(hash);
            notify_history_new(dnode->node_id,
                               dnode->sequence,
                               dnode->meta());
            ifmetrics(metrics->n_dnodes.Increment());
            pdebug("BUILD_DNODE complete: " << *dnode);
        } else {
            pdebug("BUILD_DNODE: d-node for " << hash << " already exists.");
            dnode->incr_count();
        }
        return dnode;
    }

    UnitigNode * build_unode(const std::string& sequence,
                             std::vector<hash_t>& tags,
                             hash_t left_end,
                             hash_t right_end) {

        auto lock = lock_nodes();
        id_t id = _unitig_id_counter;
        
        // Transfer the UnitigNode's ownership to the map;
        // get its new memory address
        unitig_nodes.emplace(id, std::move(make_unique<UnitigNode>(id,
                                                                   left_end,
                                                                   right_end,
                                                                   sequence)));
        UnitigNode * unode_ptr = unitig_nodes[id].get();

        _unitig_id_counter++;
        _n_unitig_nodes++;
        _n_updates++;
        ifmetrics(metrics->n_unodes.Increment());

        // Link up its new tags
        unode_ptr->tags.insert(std::end(unode_ptr->tags),
                               std::begin(tags),
                               std::end(tags));
        for (auto tag: tags) {
            unitig_tag_map.insert(make_pair(tag, unode_ptr));
        }
        unitig_end_map.insert(make_pair(left_end, unode_ptr));
        unitig_end_map.insert(make_pair(right_end, unode_ptr));

        auto unode_meta = recompute_node_meta(unode_ptr);
        unode_ptr->set_node_meta(unode_meta);
        ifmetrics(metrics->increment_cdbg_node(unode_meta));

        notify_history_new(id, unode_ptr->sequence, unode_ptr->meta());
        pdebug("BUILD_UNODE complete: " << *unode_ptr);

        return unode_ptr;
    }

    void clip_unode(direction_t clip_from,
                    hash_t old_unode_end,
                    hash_t new_unode_end) {
        
        auto lock = lock_nodes();

        auto unode = switch_unode_ends(old_unode_end, new_unode_end);
        assert(unode != nullptr);
        pdebug("CLIP: " << *unode << " from " << (clip_from == DIR_LEFT ? std::string("LEFT") : std::string("RIGHT")) <<
               " and swap " << old_unode_end << " to " << new_unode_end);

        if (unode->sequence.length() == this->_K) {
            ifmetrics(metrics->decrement_cdbg_node(unode->meta()));
            delete_unode(unode);
            pdebug("CLIP complete: deleted null unode.");
        } else {
            ifmetrics(metrics->n_clips.Increment());
            if (clip_from == DIR_LEFT) {
                unode->sequence = unode->sequence.substr(1);
                unode->set_left_end(new_unode_end);

                ifmetrics(metrics->decrement_cdbg_node(unode->meta()));
                auto meta = recompute_node_meta(unode);
                ifmetrics(metrics->increment_cdbg_node(meta));
                unode->set_node_meta(meta);

                notify_history_clip(unode->node_id, unode->sequence, unode->meta());
                pdebug("CLIP complete: " << *unode);
            } else {
                unode->sequence = unode->sequence.substr(0, unode->sequence.length() - 1);
                unode->set_right_end(new_unode_end);

                ifmetrics(metrics->decrement_cdbg_node(unode->meta()));
                auto meta = recompute_node_meta(unode);
                ifmetrics(metrics->increment_cdbg_node(meta));
                unode->set_node_meta(meta);

                notify_history_clip(unode->node_id, unode->sequence, unode->meta());
                pdebug("CLIP complete: " << *unode);
            }
        }

        ++_n_updates;
    }

    void extend_unode(direction_t ext_dir,
                      const std::string& new_sequence,
                      hash_t old_unode_end,
                      hash_t new_unode_end,
                      std::vector<hash_t>& new_tags) {

        auto lock = lock_nodes();

        auto unode = switch_unode_ends(old_unode_end, new_unode_end);
        if (unode->meta() == TRIVIAL) {
            unitig_end_map.insert(make_pair(old_unode_end, unode));
        }
        auto id = unode->node_id;

        assert(unode != nullptr); 

        pdebug("EXTEND: from " << old_unode_end << " to " << new_unode_end
               << (ext_dir == DIR_LEFT ? std::string(" to LEFT") : std::string(" to RIGHT"))
               << " adding " << new_sequence << " to"
               << std::endl << *unode);
#ifdef DEBUG_CDBG
        //auto counts = dbg->get_counts(unode->sequence);
        //for (auto c : counts) {
        //    assert(c > 0);
        //}
#endif
        
        if (ext_dir == DIR_RIGHT) {
            unode->extend_right(new_unode_end, new_sequence);
        } else {
            unode->extend_left(new_unode_end, new_sequence);
        }

        std::copy(new_tags.begin(), new_tags.end(), std::back_inserter(unode->tags));
        for (auto tag: new_tags) {
            unitig_tag_map.insert(make_pair(tag, unode));
        }

        ifmetrics(metrics->n_extends.Increment());
        ifmetrics(metrics->decrement_cdbg_node(unode->meta()));
        auto meta = recompute_node_meta(unode);
        ifmetrics(metrics->increment_cdbg_node(meta));
        unode->set_node_meta(meta);
        ++_n_updates;

        notify_history_extend(unode->node_id, unode->sequence, unode->meta());
        pdebug("EXTEND complete: " << *unode);
    }

    void split_unode(id_t node_id,
                     size_t split_at,
                     std::string split_kmer,
                     hash_t new_right_end,
                     hash_t new_left_end) {

        UnitigNode * unode;
        std::string right_unitig;
        hash_t right_unode_right_end;

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
                                  unode->sequence.substr((this->_K - 1), split_at);
                switch_unode_ends(unode->left_end(), new_left_end);
                unitig_end_map.insert(make_pair(new_right_end, unode));

                unode->set_left_end(new_left_end);
                unode->set_right_end(new_right_end);

                ifmetrics(metrics->n_splits.Increment());
                unode->set_node_meta(FULL);
                ifmetrics(metrics->decrement_cdbg_node(CIRCULAR));
                ifmetrics(metrics->increment_cdbg_node(FULL));
                ++_n_updates;

                notify_history_split_circular(unode->node_id, unode->sequence, unode->meta());
                pdebug("SPLIT complete (CIRCULAR): " << *unode);
                return;

            }
            pdebug("SPLIT: " << new_right_end << " left of root, "
                    << new_left_end << " right of root, at " << split_at
                    << std::endl << *unode);

            assert((split_at != 0) && (split_at != unode->sequence.size() - this->_K));
            right_unitig = unode->sequence.substr(split_at + 1);

            // set the left unode right end to the new right end
            right_unode_right_end = unode->right_end();
            switch_unode_ends(unode->right_end(), new_right_end);
            unode->set_right_end(new_right_end);
            unode->sequence = unode->sequence.substr(0, split_at + this->_K - 1);
            
            ifmetrics(metrics->n_splits.Increment());
            ifmetrics(metrics->decrement_cdbg_node(unode->meta()));
            auto meta = recompute_node_meta(unode);
            ifmetrics(metrics->increment_cdbg_node(meta));
            unode->set_node_meta(meta);
            ++_n_updates;
        }

        std::vector<hash_t> tags; // TODO: WARNING: broken
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
                      hash_t left_end,
                      hash_t right_end,
                      std::vector<hash_t>& new_tags) {
        /* span_sequence is the (K * 2) - 2 sequence connecting the two unitigs
         *
         */

        UnitigNode *left_unode, *right_unode;
        std::string right_sequence;
        hash_t new_right_end;

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
            ifmetrics(metrics->decrement_cdbg_node(right_unode->meta()));
            ifmetrics(metrics->n_circular_merges.Increment());
            std::string extend = span_sequence;
            //if (n_span_kmers < this->_K - 1) {
                pdebug("Overlap between merged sequence, trimming right " << n_span_kmers);
                extend = span_sequence.substr(this->_K-1, n_span_kmers);
            //} 
            extend_unode(DIR_RIGHT,
                         extend,
                         left_end, // this is left_unode's right_end
                         left_unode->left_end(),
                         new_tags);
        } else {

            pdebug("MERGE: " << left_end << " to " << right_end
                   << " with " << span_sequence 
                   << std::endl << *left_unode << std::endl << *right_unode);

            if (n_span_kmers < this->_K - 1) {
                size_t trim = (this->_K - 1) - n_span_kmers;
                pdebug("Overlap between merged sequence, trimming right " << trim);
                right_sequence = right_unode->sequence.substr(trim);
            } else {
                pdebug("No overlap, adding segment sequence, " << n_span_kmers);
                right_sequence = span_sequence.substr(this->_K - 1, n_span_kmers - this->_K + 1)
                                                       + right_unode->sequence;
            }
            std::copy(right_unode->tags.begin(), right_unode->tags.end(),
                      std::back_inserter(new_tags));
            new_right_end = right_unode->right_end();

            delete_unode(right_unode);
            extend_unode(DIR_RIGHT,
                         right_sequence,
                         left_end,
                         new_right_end,
                         new_tags);
            ifmetrics(metrics->n_merges.Increment());

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
            ifmetrics(metrics->decrement_cdbg_node(unode->meta()));
            for (hash_t tag: unode->tags) {
                unitig_tag_map.erase(tag);
            }
            unitig_end_map.erase(unode->left_end());
            unitig_end_map.erase(unode->right_end());

            unitig_nodes.erase(id);
            unode = nullptr;
            _n_unitig_nodes--;
            _n_updates++;
            ifmetrics(metrics->n_unodes.Decrement());
            ifmetrics(metrics->n_deletes.Increment());
        }
    }

    void delete_unodes_from_tags(std::vector<hash_t>& tags) {
        for (auto tag: tags) {
            UnitigNode * unode = query_unode_tag(tag);
            if (unode != nullptr) {
                delete_unode(unode);
            }
        }
    }

    void notify_history_new(id_t id, std::string& sequence, node_meta_t meta) {
        auto event = make_shared<events::HistoryNewEvent>();
        event->id = id;
        event->sequence = sequence;
        event->meta = meta;
        this->notify(event);
    }

    void notify_history_merge(id_t lparent, id_t rparent, id_t child,
                          std::string& sequence, node_meta_t meta) {
        auto event = make_shared<events::HistoryMergeEvent>();
        event->lparent = lparent;
        event->rparent = rparent;
        event->child = child;
        event->meta = meta;
        event->sequence = sequence;
        this->notify(event);
    }

    void notify_history_extend(id_t id, std::string& sequence, node_meta_t meta) {
        auto event = make_shared<events::HistoryExtendEvent>();
        event->id = id;
        event->sequence = sequence;
        event->meta = meta;
        this->notify(event);
    }

    void notify_history_clip(id_t id, std::string& sequence, node_meta_t meta) {
        auto event = make_shared<events::HistoryClipEvent>();
        event->id = id;
        event->sequence = sequence;
        event->meta = meta;
        this->notify(event);
    }

    void notify_history_split(id_t parent, id_t lchild, id_t rchild,
                          std::string& lsequence, std::string& rsequence,
                          node_meta_t lmeta, node_meta_t rmeta) {
        auto event = make_shared<events::HistorySplitEvent>();
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
        auto event = make_shared<events::HistorySplitCircularEvent>();
        event->id = id;
        event->sequence = sequence;
        event->meta = meta;
        this->notify(event);
    }

    void validate(const std::string& filename) {
        std::ofstream out;
        out.open(filename);

        auto lock = lock_nodes();

        for (auto it = unitig_nodes.begin(); it != unitig_nodes.end(); ++it) {
            auto unode = it->second.get();
            auto counts = dbg->get_counts(unode->sequence);
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

    void write_fasta(const std::string& filename) {
        std::ofstream out;
        out.open(filename);
        write_fasta(out);
        out.close();
    }

    void write_fasta(std::ofstream& out) {
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

    void write_gfa1(const std::string& filename) {
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
                l.cigar = std::to_string(_K) + "M";

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
                l.cigar = std::to_string(_K) + "M";

                std::string link_name = get_link_name(it->second->node_id, out_node->node_id);
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

/*
template <class GraphType>
class AsyncCDBG : public cDBG<GraphType>,
                  public EventListener {

public:

    using BaseType = cDBG<GraphType>;

    AsyncCDBG(uint16_t K)
        : BaseType(K),
          EventListener("cDBG")
    {
        // Event types that the cDBG EventListener should filter for.
        // Other event types will be ignored by the EventListener
        // superclass on notify().
        //
        this->msg_type_whitelist.insert(boink::event_types::MSG_ADD_DNODE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_ADD_UNODE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_SPLIT_UNODE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_EXTEND_UNODE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_MERGE_UNODES);
        this->msg_type_whitelist.insert(boink::event_types::MSG_DELETE_UNODE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_INCR_DNODE_COUNT);
    }

    virtual void handle_msg(shared_ptr<Event> event) {
        // Superclass override for handling the filtered messages.
        // For now, unrecognized messages are a no-op.
        // Follows the general pattern for the rest of boink:
        //  static_cast the Event to its subclassed type, lock
        //  underlying datastructures as appropriate, and handle
        //  updates.
        //
        switch(event->msg_type) {
            case boink::event_types::MSG_ADD_DNODE:
                {
                    auto * data = static_cast<BuildDNodeEvent*>(event.get());
                    auto lock = this->lock_dnodes();
                    this->build_dnode(data->hash, data->kmer);
                }
                return;
            case boink::event_types::MSG_ADD_UNODE:
                {
                    auto * data = static_cast<BuildUNodeEvent*>(event.get());
                    auto lock = this->lock_unodes();
                    this->build_unode(data->sequence,
                                      data->tags,
                                      data->left_end,
                                      data->right_end);
                }
                return;
            case boink::event_types::MSG_DELETE_UNODE:
                {
                    auto * data = static_cast<DeleteUNodeEvent*>(event.get());
                    auto lock = this->lock_unodes();
                    this->delete_unode(this->query_unode_id(data->node_id));
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
        auto lock = this->lock_dnodes();
        if (this->query_dnode(hash) == nullptr) {
            this->decision_nodes.insert(make_pair(hash, nullptr));
        }
        pdebug("Built d-node marker: " << hash);
    }

    bool query_dnode_marker(hash_t hash) {
        // sync query for d-nodes that are either
        // constructed already or marked for construction
        // synchronously
        auto search = this->decision_nodes.find(hash);
        if (search != this->decision_nodes.end()) {
            return true;
        }
        return false;       
    }
};
*/
}
}

#undef pdebug
#endif
