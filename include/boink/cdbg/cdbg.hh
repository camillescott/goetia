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
#include <gfakluge.hpp>
#include "boink/storage/sparsepp/spp.h"
#pragma GCC diagnostic pop

#include "boink/boink.hh"

#include "boink/events.hh"
#include "boink/event_types.hh"
#include "boink/metrics.hh"
#include "boink/reporting/reporters.hh"

#include "boink/assembly.hh"
#include "boink/hashing/hashing_types.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/kmers/kmerclient.hh"
#include "boink/minimizers.hh"
#include "boink/storage/storage.hh"

#include <prometheus/registry.h>

#include "boink/cdbg/cdbg_types.hh"
#include "boink/cdbg/metrics.hh"\

#include "utils/stringutils.h"

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


template <class GraphType>
struct cDBG {


    protected:

        using ShifterType   = typename GraphType::shifter_type;
        using TraversalType = Traverse<GraphType>;
        using MinimizerType = typename WKMinimizer<ShifterType>::Minimizer;

    public:

    typedef GraphType                         graph_type;

    typedef ShifterType                       shifter_type;
    typedef typename shifter_type::hash_type  hash_type;
    typedef typename shifter_type::kmer_type  kmer_type;
    typedef typename shifter_type::shift_type shift_type;

    typedef TraversalType traverser_type;
    typedef MinimizerType minimizer_type;


    class Graph : public kmers::KmerClient,
                  public events::EventNotifier {

        /* Map of k-mer hash --> DecisionNode. DecisionNodes take
         * their k-mer hash value as their Node ID.
         */
        typedef spp::sparse_hash_map<hash_type,
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

        // The actual k-mer hash --> DNode map
        dnode_map_t decision_nodes;

        // The actual ID --> UNode map
        unode_map_t unitig_nodes;
        // The map from Unitig end k-mer hashes to UnitigNodes
        spp::sparse_hash_map<hash_type, UnitigNode*> unitig_end_map;
        // The map from dBG k-mer tags to UnitigNodes
        spp::sparse_hash_map<hash_type, UnitigNode*> unitig_tag_map;

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

        std::shared_ptr<prometheus::Registry> pr_registry;

    public:

        std::shared_ptr<GraphType> dbg;
        std::shared_ptr<cDBGMetrics> metrics;

        Graph(std::shared_ptr<GraphType> dbg,
              std::shared_ptr<prometheus::Registry> metrics_registry,
              uint64_t minimizer_window_size=8);

        std::unique_lock<std::mutex> lock_nodes() {
            return std::unique_lock<std::mutex>(mutex);
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

        uint64_t n_unitig_ends() const {
            return unitig_end_map.size();
        }

        /* Node query methods: separate query mechanisms for
         * decision nodes and unitig nodes.
         */

        CompactNode* query_cnode(hash_type hash);

        DecisionNode* query_dnode(hash_type hash);

        vector<DecisionNode*> query_dnodes(const std::string& sequence);

        UnitigNode * query_unode_end(hash_type end_kmer);

        UnitigNode * query_unode_tag(hash_type hash);

        UnitigNode * query_unode_id(id_t id);

        bool has_dnode(hash_type hash);

        bool has_unode_end(hash_type end_kmer);

        CompactNode * find_rc_cnode(CompactNode * root);

        /* Neighbor-finding and traversal.
         *
         */

        std::pair<std::vector<CompactNode*>,
                  std::vector<CompactNode*>> find_dnode_neighbors(DecisionNode* dnode);

        std::pair<DecisionNode*, DecisionNode*> find_unode_neighbors(UnitigNode * unode);

        vector<CompactNode*> traverse_breadth_first(CompactNode* root);

        spp::sparse_hash_map<id_t, std::vector<id_t>> find_connected_components();

        /*
         * Graph Mutation
         */

        node_meta_t recompute_node_meta(UnitigNode * unode);

        UnitigNode * switch_unode_ends(hash_type old_unode_end,
                                       hash_type new_unode_end);

        DecisionNode* build_dnode(hash_type hash,
                                  const std::string& kmer);

        UnitigNode * build_unode(const std::string& sequence,
                                 std::vector<hash_type>& tags,
                                 hash_type left_end,
                                 hash_type right_end);

        void clip_unode(direction_t clip_from,
                        hash_type old_unode_end,
                        hash_type new_unode_end);

        void extend_unode(direction_t ext_dir,
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

        void delete_unode(UnitigNode * unode);
    
        void delete_unodes_from_tags(std::vector<hash_type>& tags);

        void delete_dnode(DecisionNode * dnode);

        /*
         * Event notification
         */

        void notify_history_new(id_t id, std::string& sequence, node_meta_t meta);

        void notify_history_merge(id_t lparent, id_t rparent, id_t child,
                                  std::string& sequence, node_meta_t meta);

        void notify_history_extend(id_t id, std::string& sequence, node_meta_t meta);

        void notify_history_clip(id_t id, std::string& sequence, node_meta_t meta);

        void notify_history_split(id_t parent, id_t lchild, id_t rchild,
                                  std::string& lsequence, std::string& rsequence,
                                  node_meta_t lmeta, node_meta_t rmeta);

        void notify_history_split_circular(id_t id, std::string& sequence, node_meta_t meta);


        /*
         * File output
         */

        void validate(const std::string& filename);

        void write(const std::string& filename, cDBGFormat format);

        void write(std::ofstream& out, cDBGFormat format);

        void write_fasta(const std::string& filename);

        void write_fasta(std::ofstream& out);

        void write_gfa1(const std::string& filename);

        void write_gfa1(std::ofstream& out);

        void write_graphml(const std::string& filename,
                           const std::string graph_name="cDBG");

        void write_graphml(std::ofstream& out,
                           const std::string graph_name="cDBG");
    };


    class ComponentReporter : public boink::reporting::SingleFileReporter {
    public:

        class Metrics {

            private:

                prometheus::Family<prometheus::Summary>&   recompute_time_family;
                prometheus::Summary::Quantiles             recompute_time_quantiles;
                
            public:
                using Quantile = prometheus::detail::CKMSQuantiles::Quantile;

                prometheus::Summary&                       recompute_time;

            private:

                prometheus::Family<prometheus::Gauge>&     component_counts_family;

            public:

                prometheus::Gauge&                         n_components;
                prometheus::Gauge&                         max_component_size;
                prometheus::Gauge&                         min_component_size;

                Metrics(std::shared_ptr<prometheus::Registry> registry)
                    : recompute_time_family     (prometheus::BuildSummary()
                                                             .Name("boink_cdbg_components_compute_time_seconds")
                                                             .Register(*registry)),
                      recompute_time_quantiles  {{.75, .1},
                                                 {.5,  .1},
                                                 {.25, .1}},
                      recompute_time            (recompute_time_family.Add({{"time", "quantiles"}},
                                                                           recompute_time_quantiles)),
                      component_counts_family   (prometheus::BuildGauge()
                                                             .Name("boink_cdbg_components_current_total")
                                                             .Register(*registry)),
                      n_components              (component_counts_family.Add({{"size", "all_components"}})),
                      max_component_size        (component_counts_family.Add({{"size", "max_component"}})),
                      min_component_size        (component_counts_family.Add({{"size", "min_component"}}))
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
                          std::shared_ptr<prometheus::Registry>  registry,
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

            metrics = std::make_unique<ComponentReporter::Metrics>(registry);
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

            metrics->n_components.Set(components.size());
            metrics->max_component_size.Set(max_component);
            metrics->min_component_size.Set(min_component);

            auto time_elapsed = std::chrono::system_clock::now() - time_start;
            metrics->recompute_time.Observe(std::chrono::duration<double>(time_elapsed).count());
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

        virtual void handle_msg(std::shared_ptr<events::Event> event) {
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
                   const string& output_prefix)
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
                                             const string& output_prefix) {
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
}

#undef pdebug
#endif
