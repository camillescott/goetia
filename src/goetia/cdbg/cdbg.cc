/**
 * (c) Camille Scott, 2019
 * File   : cdbg.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 03.09.2019
 */

#include "goetia/cdbg/cdbg.hh"

#include "goetia/dbg.hh"
#include "goetia/hashing/rollinghashshifter.hh"
#include "goetia/storage/storage_types.hh"

// save diagnostic state
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wchar-subscripts"
#include "goetia/parsing/gfakluge/gfakluge.hpp"
#pragma GCC diagnostic pop


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


template <template <class, class> class GraphType,
          class StorageType,
          class ShifterType>
cDBG<GraphType<StorageType, ShifterType>>::
CompactNode::CompactNode(id_t node_id,
                    const std::string& sequence,
                    node_meta_t meta)
            : _meta(meta),
              node_id(node_id),
              component_id(NULL_ID),
              sequence(sequence)
        {
        }


template <template <class, class> class GraphType,
          class StorageType,
          class ShifterType>
cDBG<GraphType<StorageType, ShifterType>>::
DecisionNode::DecisionNode(id_t node_id, const std::string& sequence)
            : CompactNode(node_id, sequence, DECISION),
              _dirty(true),
              _left_degree(0),
              _right_degree(0),
              _count(1)
        {    
        }


template <template <class, class> class GraphType,
          class StorageType,
          class ShifterType>
cDBG<GraphType<StorageType, ShifterType>>::
UnitigNode::UnitigNode(id_t node_id,
                   hash_type left_end,
                   hash_type right_end,
                   const std::string& sequence,
                   node_meta_t meta)
            : CompactNode(node_id, sequence, meta),
              _left_end(left_end),
              _right_end(right_end) { 
        }

template <template <class, class> class GraphType,
          class StorageType,
          class ShifterType>
cDBG<GraphType<StorageType, ShifterType>>::
Graph::Graph(std::shared_ptr<graph_type> dbg,
             uint64_t minimizer_window_size)
    : K(dbg->K),
      dbg(dbg),
      _n_updates(0),
      _unitig_id_counter(UNITIG_START_ID),
      _n_unitig_nodes(0),
      component_id_counter(0)
{
    metrics = std::make_shared<cDBGMetrics>();
}


template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
auto
cDBG<GraphType<StorageType, ShifterType>>::
Graph::find_unode_neighbors(UnitigNode * unode)
-> std::pair<DecisionNode*, DecisionNode*> {

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


template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
auto
cDBG<GraphType<StorageType, ShifterType>>::
Graph::find_dnode_neighbors(DecisionNode* dnode)
-> std::pair<std::vector<CompactNode*>,
             std::vector<CompactNode*>> {

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


template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
auto
cDBG<GraphType<StorageType, ShifterType>>::
Graph::traverse_breadth_first(CompactNode* root)
-> std::vector<CompactNode*> {

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


template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
auto
cDBG<GraphType<StorageType, ShifterType>>::
Graph::find_connected_components()
-> spp::sparse_hash_map<id_t, std::vector<id_t>>{

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


template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
auto
cDBG<GraphType<StorageType, ShifterType>>::
Graph::recompute_node_meta(UnitigNode * unode)
-> node_meta_t {

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



/*
 * Graph Mutation methods
 *
 */

template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
auto
cDBG<GraphType<StorageType, ShifterType>>::
Graph::switch_unode_ends(hash_type old_unode_end,
                         hash_type new_unode_end)
-> UnitigNode * {

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


template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
auto
cDBG<GraphType<StorageType, ShifterType>>::
Graph::build_dnode(hash_type hash,
                   const std::string& kmer)
-> DecisionNode* {

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
        metrics->n_dnodes++;
        pdebug("BUILD_DNODE complete: " << *dnode);
    } else {
        pdebug("BUILD_DNODE: d-node for " << hash << " already exists.");
        dnode->incr_count();
    }
    return dnode;
}


template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
auto
cDBG<GraphType<StorageType, ShifterType>>::
Graph::build_unode(const std::string& sequence,
                   std::vector<hash_type>& tags,
                   hash_type left_end,
                   hash_type right_end)
-> UnitigNode * {

    auto lock = lock_nodes();
    id_t id = _unitig_id_counter;

    if (ShifterType::hash(sequence.substr(0, this->K), this->K) != left_end) {
        std::cout << "left: " << left_end << std::endl;
        std::cout << "right: " << right_end << std::endl;
        std::cout << "seq: " << sequence << std::endl;
        std::cout << "hashed: " << ShifterType::hash(sequence.substr(this->K), this->K) << std::endl;
        assert(ShifterType::hash(sequence.substr(0, this->K), this->K) == left_end);
    }
    
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

    pdebug("BUILD_UNODE complete: " << *unode_ptr);

    return unode_ptr;
}


template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
void
cDBG<GraphType<StorageType, ShifterType>>::
Graph::clip_unode(bool      clip_from,
                  hash_type old_unode_end,
                  hash_type new_unode_end) {

    auto lock = lock_nodes();

    auto unode = switch_unode_ends(old_unode_end, new_unode_end);
    assert(unode != nullptr);
    pdebug("CLIP: " << *unode << " from " << (clip_from == DIR_LEFT ? std::string("LEFT") : std::string("RIGHT")) <<
           " and swap " << old_unode_end << " to " << new_unode_end);

    if (unode->sequence.length() == this->K) {
        metrics->decrement_cdbg_node(unode->meta());
        delete_unode(unode);
        pdebug("CLIP complete: deleted null unode.");
    } else {
        metrics->n_clips++;
        if (clip_from == DIR_LEFT) {
            unode->sequence = unode->sequence.substr(1);
            unode->set_left_end(new_unode_end);

            metrics->decrement_cdbg_node(unode->meta());
            auto meta = recompute_node_meta(unode);
            metrics->increment_cdbg_node(meta);
            unode->set_node_meta(meta);

            pdebug("CLIP complete: " << *unode);
        } else {
            unode->sequence = unode->sequence.substr(0, unode->sequence.length() - 1);
            unode->set_right_end(new_unode_end);

            metrics->decrement_cdbg_node(unode->meta());
            auto meta = recompute_node_meta(unode);
            metrics->increment_cdbg_node(meta);
            unode->set_node_meta(meta);

            pdebug("CLIP complete: " << *unode);
        }
    }

    ++_n_updates;
}


template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
void
cDBG<GraphType<StorageType, ShifterType>>::
Graph::extend_unode(bool               ext_dir,
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
           << (ext_dir == DIR_LEFT ? std::string(" to LEFT") : std::string(" to RIGHT"))
           << " adding " << new_sequence << " to"
           << std::endl << *unode);
    
    if (ext_dir == DIR_RIGHT) {
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

    pdebug("EXTEND complete: " << *unode);
}



template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
void
cDBG<GraphType<StorageType, ShifterType>>::
Graph::split_unode(id_t node_id,
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

    pdebug("SPLIT complete: " << std::endl << *unode << std::endl << *new_node);

}



template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
void
cDBG<GraphType<StorageType, ShifterType>>::
Graph::merge_unodes(const std::string& span_sequence,
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
        extend_unode(DIR_RIGHT,
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
        extend_unode(DIR_RIGHT,
                     right_sequence,
                     left_end,
                     new_right_end,
                     new_tags);
        metrics->n_merges++;

    }
    
    pdebug("MERGE complete: " << *left_unode);
}



template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
void
cDBG<GraphType<StorageType, ShifterType>>::
Graph::write_gfa1(std::ofstream& out) {
    
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



/*
 * Model-level functions
 */


template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
auto
cDBG<GraphType<StorageType, ShifterType>>
::compute_connected_component_metrics(std::shared_ptr<cDBG<GraphType<StorageType, ShifterType>>::Graph> cdbg,
                                      size_t sample_size)
-> std::tuple<size_t, size_t, size_t, std::vector<size_t>> {

    auto time_start = std::chrono::system_clock::now();

    ReservoirSample<size_t> component_size_sample(sample_size);
    size_t max_component = 0;
    size_t min_component = std::numeric_limits<size_t>::max();
    auto components = cdbg->find_connected_components();

    for (auto id_comp_pair : components) {
        size_t component_size = id_comp_pair.second.size();
        component_size_sample.sample(component_size);
        max_component = (component_size > max_component) ? component_size : max_component;
        min_component = (component_size < min_component) ? component_size : min_component;
    }

    auto time_elapsed = std::chrono::system_clock::now() - time_start;
    _cerr("Finished recomputing components. Elapsed time: " <<
          std::chrono::duration<double>(time_elapsed).count());

    return {components.size(), min_component, max_component, component_size_sample.get_result()};
}


template <template <class, class> class GraphType,
          class StorageType, 
          class ShifterType>
auto
cDBG<GraphType<StorageType, ShifterType>>
::compute_unitig_fragmentation(std::shared_ptr<Graph> cdbg,
                               std::vector<size_t>    bins)
-> std::vector<size_t> {

    auto time_start = std::chrono::system_clock::now();
    auto lock       = cdbg->lock_nodes();
    _cerr("Summing unitig length bins...");

    std::vector<size_t> bin_sums(bins.size(), 0);

    for (auto it = cdbg->unodes_begin(); it != cdbg->unodes_end(); ++it) {
        auto seq_len = it->second->sequence.length();
        for (size_t bin_num = 0; bin_num < bins.size() - 1; bin_num++) {
            if (seq_len >= bins[bin_num] && seq_len < bins[bin_num+1]) {
                bin_sums[bin_num] += (seq_len - cdbg->K + 1);
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

template class cDBG<goetia::dBG<BitStorage, FwdLemireShifter>>;
template class cDBG<goetia::dBG<BitStorage, CanLemireShifter>>;

template class cDBG<goetia::dBG<SparseppSetStorage, FwdLemireShifter>>;
template class cDBG<goetia::dBG<SparseppSetStorage, CanLemireShifter>>;

template class cDBG<goetia::dBG<PHMapStorage, FwdLemireShifter>>;
template class cDBG<goetia::dBG<PHMapStorage, CanLemireShifter>>;

template class cDBG<goetia::dBG<BTreeStorage, FwdLemireShifter>>;
template class cDBG<goetia::dBG<BTreeStorage, CanLemireShifter>>;

template class cDBG<goetia::dBG<ByteStorage, FwdLemireShifter>>;
template class cDBG<goetia::dBG<ByteStorage, CanLemireShifter>>;

template class cDBG<goetia::dBG<NibbleStorage, FwdLemireShifter>>;
template class cDBG<goetia::dBG<NibbleStorage, CanLemireShifter>>;

template class cDBG<goetia::dBG<QFStorage, FwdLemireShifter>>;
template class cDBG<goetia::dBG<QFStorage, CanLemireShifter>>;

}
