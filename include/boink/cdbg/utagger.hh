/**
 * (c) Camille Scott, 2019
 * File   : utagger.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 05.08.2019
 */

#ifndef BOINK_UTAGGER_HH
#define BOINK_UTAGGER_HH

#include <memory>
#include <optional>

#include "boink/processors.hh"
#include "boink/traversal.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/hashing/hash_combine.hh"
#include "boink/dbg.hh"
#include "boink/cdbg/udbg.hh"
#include "boink/storage/storage_types.hh"

namespace boink {
namespace cdbg {

template <class StorageType>
struct UTagger {

    typedef dBG<StorageType, hashing::CanUnikmerShifter> graph_type;
    typedef dBGWalker<graph_type>                        walker_type;

    // inject dependent typename boilerplate: see boink/meta.hh
    _boink_model_typedefs_from_graphtype(graph_type);
    _boink_walker_typedefs_from_graphtype(walker_type);

    typedef typename shifter_type::ukhs_type   ukhs_type;
    typedef std::tuple<value_type, value_type> link_type;

    struct Tag {
        uint32_t left_partition;
        uint32_t right_partition;
        //link_type  link;
    };

    class StreamingTagger :public events::EventNotifier {

    public:

        typedef spp::sparse_hash_map<link_type,
                                     Tag,
                                     hash_tuple::hash<link_type>> tag_map_t;
        typedef typename tag_map_t::const_iterator                tag_map_iter_t;

        typedef spp::sparse_hash_set<link_type,
                                     hash_tuple::hash<link_type>> link_set_t;
        typedef typename link_set_t::const_iterator               link_set_iter_t;

        tag_map_t                                                 tag_map;

        std::shared_ptr<graph_type> dbg;
        std::shared_ptr<ukhs_type>  ukhs;

        const uint16_t              W;
        const uint16_t              K;

        StreamingTagger(std::shared_ptr<graph_type>& dbg,
                        std::shared_ptr<ukhs_type>&  ukhs)
            : EventNotifier(),
              dbg(dbg),
              ukhs(ukhs),
              W(dbg->K),
              K(ukhs->K)
        {
        }

        static std::shared_ptr<StreamingTagger> build(std::shared_ptr<graph_type> dbg,
                                                      std::shared_ptr<ukhs_type>  ukhs) {
            return std::make_shared<StreamingTagger>(dbg, ukhs);
        }

        /**
         * @Synopsis  Creates all tags in sequence's neighborhood and adds them to the master 
         *            tag map, while adding the hashes from sequence to the underlying de Bruijn graph.
         *
         * @Param sequence Sequence to insert.
         *
         * @Returns   Number of tags created.
         */
        size_t insert_sequence(const std::string& sequence) {
            auto new_tags = find_new_tags(sequence);
            if (new_tags.size()) {
                for (const auto& item : new_tags) {
                    tag_map[item.first] = item.second;
                }
                return new_tags.size();
            }
            
            return 0;
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
            hashing::KmerIterator<graph_type> kmer_iter(sequence, dbg.get());
            std::vector<std::tuple<hash_type, size_t, shift_pair_type>> results;

            size_t position = 0;
            while(!kmer_iter.done()) {
                auto hash = kmer_iter.next();
                if (dbg->insert(hash)) {
                    // Note that left_extensions() and right_extensions() only collect
                    // up the *potential* neighbors for the cursor. We do this here
                    // for all new k-mers in order to minimize the amount of hashing
                    // we have to do; actually querying them against the graph and filtering
                    // them out has to wait until all the new k-mers from the sequence
                    // have been inserted.
                    results.emplace_back(hash, position, std::make_pair(dbg->left_extensions(),
                                                                        dbg->right_extensions()));
                } 
                ++position;
            }

            return std::move(results);
        }


        /**
         * @Synopsis  Finds all new tags induced by the insertion of sequence and returns them.
         *
         * @Param sequence 
         *
         * @Returns   The new tags.
         */
        tag_map_t find_new_tags(const std::string& sequence) {

            auto hashes = find_new_extensions(sequence);

            tag_map_t tags;
            if (hashes.size() == 0) {
                return tags;
            }
 
            for (const auto& hash_tuple : hashes) {
                // Now that all the new k-mers have been inserted, we can filter the nodes
                // against the dBG and reduce the collection down to only the actual neighbors.
                create_neighborhood_tags(std::get<0>(hash_tuple),
                                         dbg->filter_nodes(std::get<2>(hash_tuple)),
                                         tags);
            }

            return tags;
        }

        /**
         * @Synopsis  Creates all possible tags in the given neighbood surrounding root
         *            and adds them to the map given by tags.
         *
         * @Param root      The neighborhood root.
         * @Param neighbors The (filtered) neighbors of root.
         * @Param tags      Map to store tags in.
         */
        void create_neighborhood_tags(const hash_type&       root,
                                      const shift_pair_type& neighbors,
                                      tag_map_t&             tags) {
            // check in-neighbors
            for (const shift_type<hashing::DIR_LEFT>& in_neighbor : neighbors.first) {
                auto tag = create_tag(in_neighbor, root);
                if (tag) {
                    tags[std::make_pair(in_neighbor, root)] = std::move(tag.value());
                }
            }
            // and out neighbors
            for (const shift_type<hashing::DIR_RIGHT>& out_neighbor : neighbors.second) {
                auto tag = create_tag(root, out_neighbor);
                if (tag) {
                    tags[std::make_pair(root, out_neighbor)] = std::move(tag.value());
                }
            }

        }

        /**
         * @Synopsis  Checks if (first, second) could be a Tag, and if so, creates the appropriate
         *            tag and returns it.
         *
         * @Param first
         * @Param second
         *
         * @Returns   An std::optional<Tag> containing the Tag if successful.
         */
        std::optional<Tag> create_tag(const hash_type& first, const hash_type& second) {
            if (is_tag(first, second)) {
                //std::cerr << "\tunikmers don't match" << std::endl;
                //std::cerr << "\tfirst: " << first << ", second:" << second << std::endl;
                Tag tag;
                tag.left_partition = first.minimizer.partition;
                tag.right_partition = second.minimizer.partition;
                //tag.link = std::make_pair(first.hash.value(),
                //                          second.hash.value());
                return tag;
            } else {
                return {};
            }
        }

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

        /**
         * @Synopsis  Check if the pair (u,v) are a tag link and return the corresponding tag if so.
         *
         * @Param u 
         * @Param v
         *
         * @Returns An std::optional<Tag> containing the tag if found.  
         */
        std::optional<Tag> query_tag(const hash_type& u, const hash_type& v) {
            if (!is_tag(u, v)) {
                return {};
            }
            return query_tag(std::make_pair(u.value(), v.value()));
        }

        std::optional<Tag> query_tag(const link_type& link) {
            auto search = tag_map.find(link);
            if (search != tag_map.end()) {
                return search->second;
            }

            return {};
        }

        /**
         * @Synopsis  Find tags within the given sequence.
         *
         * @Param sequence Sequence to search.
         *
         * @Returns   Ordered vector of tags from sequence.
         */
        std::vector<std::pair<link_type, Tag>> query_sequence_tags(const std::string& sequence) {

            std::vector<std::pair<link_type, Tag>> tags;
            hashing::KmerIterator<graph_type> kmer_iter(sequence, dbg.get());

            hash_type u = kmer_iter.next();
            if (kmer_iter.done()) {
                return tags;
            }
            while(!kmer_iter.done()) {
                hash_type v = kmer_iter.next();
                if (is_tag(u, v)) {
                    auto link = std::make_pair(u.value(), v.value());
                    auto tag = query_tag(link);
                    if (tag) {
                        tags.emplace_back(link, tag.value());
                    }
                }
                u = v;
            }

            return std::move(tags);
        }
                            
        /**
         * @Synopsis  Find (unique) tags within the given sequence or directly adjacent to its k-mers.
         *
         * @Param sequence Sequence to search.
         *
         * @Returns   Set of links corresponding to neighboring tags.
         */
        tag_map_t query_neighborhood_tags(const std::string& sequence) {

            hashing::KmerIterator<graph_type> kmer_iter(sequence, dbg.get());
            tag_map_t found;

            while(!kmer_iter.done()) {
                hash_type root = kmer_iter.next();
                auto neighborhood = dbg->filter_nodes(std::make_pair(dbg->left_extensions(),
                                                                     dbg->right_extensions()));
                for (const shift_type<hashing::DIR_LEFT>& in_neighbor : neighborhood.first) {
                    if (auto tag = query_tag(in_neighbor, root)) {
                        found[std::make_pair(in_neighbor, root)] = tag.value();
                    }
                }
                // and out neighbors
                for (const shift_type<hashing::DIR_RIGHT>& out_neighbor : neighborhood.second) {
                    if (auto tag = query_tag(root, out_neighbor)) {
                        found[std::make_pair(root, out_neighbor)] = tag.value();
                    }
                }
            }

            return std::move(found);
        }

    };


    using Processor = InserterProcessor<StreamingTagger>;

};

extern template class cdbg::UTagger<storage::BitStorage>;
extern template class cdbg::UTagger<storage::ByteStorage>;
extern template class cdbg::UTagger<storage::NibbleStorage>;
extern template class cdbg::UTagger<storage::QFStorage>;
extern template class cdbg::UTagger<storage::SparseppSetStorage>;

}
}

#endif
