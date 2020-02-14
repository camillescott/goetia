/**
 * (c) Camille Scott, 2020
 * File   : usparsetagger.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 11.02.2020
 */

#ifndef BOINK_USPARSETAGGER_HH
#define BOINK_USPARSETAGGER_HH

#include <memory>
#include <optional>

#include "boink/dbg.hh"
#include "boink/hashing/hash_combine.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/processors.hh"
#include "boink/rx_ranges.hpp"
#include "boink/storage/storage_types.hh"
#include "boink/traversal.hh"

namespace boink::cdbg {

    template <typename StorageType>
    struct USparseGraph {
        typedef dBG<StorageType, hashing::CanUnikmerShifter> graph_type;
        typedef dBGWalker<graph_type>                        walker_type;

        // inject dependent typename boilerplate: see boink/meta.hh
        _boink_model_typedefs_from_graphtype(graph_type);
        _boink_walker_typedefs_from_graphtype(walker_type);

        typedef typename shifter_type::ukhs_type   ukhs_type;
        typedef std::tuple<value_type, value_type> link_type;

        struct Tip {
            kmer_type  kmer;
            value_type partner;
            value_type neighbor_tag;
            bool       position;
        };

        struct Tag {
            value_type lneighbor, rneighbor;
        };

        typedef spp::sparse_hash_map<value_type, Tag> tag_map_t;
        typedef typename tag_map_t::const_iterator    tag_map_iter_t;

        typedef spp::sparse_hash_map<value_type, Tip> tip_map_t;
        typedef typename tip_map_t::const_iterator    tip_map_iter_t;


        class Graph : public events::EventNotifier {
          public:
            std::shared_ptr<graph_type> dbg;
            std::shared_ptr<ukhs_type>  ukhs;
            tag_map_t                   tags;
            tip_map_t                   tips;

            const uint16_t W;
            const uint16_t K;

            hash_tuple::hash<link_type> hash_link;

            struct TagStopper {
                std::optional<hash_type> prev;
                std::optional<link_type> found_link;
                Graph*                   sparse_graph;

                explicit TagStopper(Graph* sparse_graph) : prev{}, found_link{}, sparse_graph(sparse_graph) {
                }

                bool operator()(hash_type& node) {
                    if (prev && sparse_graph->is_tag(prev.value(), node)) {
                        found_link = std::make_tuple(prev.value().value(), node.value());
                        return false;
                    }
                    prev = node;
                    return true;
                }
            };

            Graph(std::shared_ptr<graph_type>& dbg, std::shared_ptr<ukhs_type>& ukhs)
                : EventNotifier(), dbg(dbg), ukhs(ukhs), W(dbg->K), K(ukhs->K) {
            }

            static std::shared_ptr<Graph> build(std::shared_ptr<graph_type>& dbg, std::shared_ptr<ukhs_type>& ukhs) {
                return std::make_shared<Graph>(dbg, ukhs);
            }

            /**
             * @brief Hashes each k-mer in sequence, and gets the set of extensions for each
             *        k-mer that is new to the dBG.
             *
             * @param sequence
             * @return std::vector<std::tuple<hash_type,
             *                                size_t,
             *                                shift_pair_type>>  Each tuple contains the hash, position in sequence,
             *                                                   and extensions, respectively.
             */
            auto find_new_extensions(const std::string& sequence)
                -> std::vector<std::tuple<hash_type, size_t, shift_pair_type>> {

                // don't confuse dbg.get(), retrieves the raw pointer,
                // with dbg->get(), which reaches through and gets the cursor hash value
                hashing::KmerIterator<graph_type>                           kmer_iter(sequence, dbg.get());
                std::vector<std::tuple<hash_type, size_t, shift_pair_type>> results;

                size_t position = 0;
                while (!kmer_iter.done()) {
                    auto hash = kmer_iter.next();
                    if (dbg->insert(hash)) {
                        // Note that left_extensions() and right_extensions() only collect
                        // up the *potential* neighbors for the cursor. We do this here
                        // for all new k-mers in order to minimize the amount of hashing
                        // we have to do; actually querying them against the graph and filtering
                        // them out has to wait until all the new k-mers from the sequence
                        // have been inserted.
                        results.emplace_back(hash,
                                             position,
                                             std::make_pair(dbg->left_extensions(), dbg->right_extensions()));
                    }
                    ++position;
                }

                return std::move(results);
            }

            /**
             * @brief Call filter_nodes on each extension in extensions to reduce them down to
             *        actual neighbors in the dBG.
             *
             * @param extensions
             * @return std::vector<std::tuple<hash_type, size_t, shift_pair_type>>  The filtered extensions.
             */
            auto filter_new_extensions(std::vector<std::tuple<hash_type, size_t, shift_pair_type>>& extensions)
                -> std::vector<std::tuple<hash_type, size_t, shift_pair_type>> {

                for (auto& ext_tuple : extensions) {
                    std::get<2>(ext_tuple) = dbg->filter_nodes(std::get<2>(ext_tuple));
                }

                return std::move(extensions);
            }

            /**
             * @brief Given the list of dbg-filtered neighbors, split them up into chains of connected
             *        new-kmers.
             *
             * @param filtered
             * @return std::vector<std::vector<std::tuple<hash_type, size_t, shift_pair_type>>>
             */
            auto build_new_segments(std::vector<std::tuple<hash_type, size_t, shift_pair_type>>& filtered)
                -> std::vector<std::vector<std::tuple<hash_type, size_t, shift_pair_type>>> {

                std::vector<std::vector<std::tuple<hash_type, size_t, shift_pair_type>>> segments;
                segments.emplace_back();  // make sure there's an empty space to add to

                if (filtered.size() == 1) {
                    segments.back().push_back(std::move(filtered.front()));

                    return std::move(segments);
                }

                size_t i = 1;
                while (i < filtered.size()) {
                    auto u = filtered.at(i - 1);
                    auto v = filtered.at(i);

                    segments.back().push_back(std::move(u));

                    if (std::get<1>(v) != std::get<1>(u) + 1) {
                        segments.emplace_back();
                    }

                    ++i;
                }

                segments.back().push_back(std::move(filtered.back()));

                return std::move(segments);
            }

            /*
                        auto build_segments(const std::string&                                           sequence,
                                            std::vector<std::tuple<hash_type, size_t, shift_pair_type>>& extensions)
                            -> std::vector<std::tuple<kmer_type,
                                                      shift_pair_type,
                                                      size_t,  //
                                                      kmer_type,
                                                      shift_pair_type,
                                                      size_t>> {

                            std::vector<std::tuple<kmer_type,
                                                   shift_pair_type,
                                                   size_t,  // left and right tips
                                                   kmer_type,
                                                   shift_pair_type,
                                                   size_t>>
                                segments;

                            std::optional<std::tuple<hash_type, size_t, shift_pair_type>> left = extensions.front();

                            size_t i = 1;
                            while (i < extensions.size()) {
                                auto u = extensions.at(i - 1);
                                auto v = extensions.at(i);

                                auto u_position  = std::get<1>(u);
                                auto u_neighbors = std::get<2>(u);
                                auto v_position  = std::get<1>(v);
                                auto v_neighbors = std::get<2>(v);

                                if (!left) {
                                    left = v;
                                }

                                if (left) {
                                    if ((neighbors.first.size() == 1 || neighbors.first.size() == 0) &&
                                        neighbors.second.size() > 1) {
                                        // right partner is a FWD decision node
                                        segments.emplace_back(std::get<0>(left.value()),
                                                              std::get<2>(left.value()),
                                                              std::get<1>(left.value()),
                                                              std::get<0>(u),
                                                              u_neighbors,
                                                              u_position);
                                    } else if (v_position != u_position + 1         // we moved to a different segment
                                               || (v_neighbors.first.size() > 1 &&  // v is a REV decision nodes
                                                   (v_neighbors.second.size() == 1 || v_neighbors.second.size() == 0)))
               { segments.emplace_back(std::get<0>
                                    }

                                    left = {};
                                }

                                // entered new segment positionally
                                else if (left && position != left_position + 1) {
                                    segments.emplace_back(
                                        left.value(), left_neighbors.value(), left_position, h, neighbors, position);
                                    left = {};
                                }
                            }
                        }
                        */

            /**
             * @Synopsis  True if the pair (u,v) could be a tag; assumes they are neighbors.
             *
             * @Param u
             * @Param v
             *
             * @Returns   True if possibly a tag.
             */
            bool is_tag(const hash_type& u, const hash_type& v) {
                return u.minimizer.partition != v.minimizer.partition;
            }
        };
    };  // namespace boink::cdbg

}  // namespace boink::cdbg
#endif