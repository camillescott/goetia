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

#include "boink/traversal.hh"
#include "boink/hashing/kmeriterator.hh"
#include "boink/hashing/ukhs.hh"
#include "boink/hashing/hash_combine.hh"
#include "boink/dbg.hh"
#include "boink/cdbg/udbg.hh"
#include "boink/storage/storage.hh"

namespace boink {
namespace cdbg {

template <class StorageType>
struct UTagger {

    typedef dBG<StorageType,
                hashing::UKHS::LazyShifter>     dbg_type;
    typedef typename dbg_type::shifter_type     shifter_type;
    typedef Traverse<dbg_type>                  traversal_type;
    typedef typename traversal_type::dBG        traverser_type;

    typedef hashing::KmerIterator<shifter_type> kmer_iter_type;
    typedef typename shifter_type::hash_type    hash_type;
    typedef typename hashing::UKHS::value_type  value_type;
    typedef typename shifter_type::kmer_type    kmer_type;
    typedef typename shifter_type::shift_type   shift_type;
    typedef TraversalState::State               state_type;

    typedef std::tuple<value_type, value_type>  link_type;

    struct Tag {
        value_type left_partition;
        value_type right_partition;
        link_type  link;
    };

    class StreamingTagger : public traverser_type,
                            public events::EventNotifier {

    public:

        using traverser_type::filter_nodes;
        using traverser_type::find_left_kmers;
        using traverser_type::find_right_kmers;
        using traverser_type::gather_left;
        using traverser_type::gather_right;
        using traverser_type::traverse_left;
        using traverser_type::traverse_right;
        using traverser_type::get_decision_neighbors;

        typedef spp::sparse_hash_map<link_type,
                                     Tag,
                                     hash_tuple::hash<link_type>> tag_map_t;
        typedef typename tag_map_t::const_iterator                tag_map_iter_t;

        tag_map_t                                                 tag_map;

        std::shared_ptr<dbg_type>            dbg;
        std::shared_ptr<hashing::UKHS::Map>  ukhs;
        const uint16_t                       unikmer_k;

        StreamingTagger(std::shared_ptr<dbg_type> dbg,
                        std::shared_ptr<hashing::UKHS::Map> ukhs)
            : traverser_type(dbg->K(), ukhs->K(), ukhs),
              EventNotifier(),
              dbg(dbg),
              ukhs(ukhs),
              unikmer_k(ukhs->K())
        {
        }

        size_t insert_sequence(const std::string& sequence) {
            std::deque<hash_type> hashes;
            auto new_tags = find_new_tags(sequence, hashes);
            for (const auto& item : new_tags) {
                tag_map[item.first] = item.second;
            }

            return new_tags.size();
        }

        tag_map_t find_new_tags(const std::string&     sequence,
                                std::deque<hash_type>& hashes) {


            kmer_iter_type kmer_iter(sequence, this);
            std::set<hash_type> new_kmers;
            std::deque<std::pair<std::vector<shift_type>,
                                 std::vector<shift_type>>> neighbors;
            while(!kmer_iter.done()) {
                auto hash = kmer_iter.next();
                hashes.push_back(hash);
                if (dbg->insert(hash)) {
                    // If the k-mer is new, track it and gather it's potential
                    // neighbors to save on hashing later
                    new_kmers.insert(hash);
                    neighbors.push_back(std::make_pair(this->gather_left(),
                                                       this->gather_right()));
                }
            }

            tag_map_t tags; 
            bool first_new = new_kmers.count(hashes.front()), second_new;
            auto first = hashes.cbegin();
            auto second = hashes.cbegin() + 1;
            while (second != hashes.cend()) {

                //std::cerr << "----------" << std::endl;
                //std::cerr << "first: " << *first << ", second: " << *second << std::endl;

                second_new = new_kmers.count(*second);
                
                if (first_new || second_new) {
                    //std::cerr << "  checking internal pair." << std::endl;
                    auto tag = create_tag(*first, *second);
                    if (tag) {
                        tags[tag.value().link] = std::move(tag.value());
                    }
                }

                if (first_new) {
                    //std::cerr << "  checking neighbors." << std::endl;
                    find_neighborhood_tags(*first,
                                           this->filter_nodes(dbg.get(),
                                                              neighbors.front()),
                                           tags);
                    neighbors.pop_front();
                }

                first_new = second_new;

                ++first;
                ++second;
            }

            if (second_new) {
                find_neighborhood_tags(hashes.back(),
                                       this->filter_nodes(dbg.get(),
                                                          neighbors.front()),
                                       tags);
            }

            return tags;
        }

        void find_neighborhood_tags(const hash_type&                          root,
                                    const std::pair<std::vector<shift_type>,
                                                    std::vector<shift_type>>& neighbors,
                                    tag_map_t&                                tags) {
            // check in-neighbors
            for (auto& in_neighbor : neighbors.first) {
                auto tag = create_tag(in_neighbor.hash, root);
                if (tag) {
                    tags[tag.value().link] = std::move(tag.value());
                }
            }
            // and out neighbors
            for (auto& out_neighbor : neighbors.second) {
                auto tag = create_tag(root, out_neighbor.hash);
                if (tag) {
                    tags[tag.value().link] = std::move(tag.value());
                }
            }

        }

        
        std::optional<Tag> create_tag(const hash_type& first, const hash_type& second) {
            if (first.unikmer.partition != second.unikmer.partition) {
                //std::cerr << "\tunikmers don't match" << std::endl;
                //std::cerr << "\tfirst: " << first << ", second:" << second << std::endl;
                Tag tag;
                tag.left_partition = first.unikmer.partition;
                tag.right_partition = second.unikmer.partition;
                tag.link = std::make_pair(first.hash,
                                          second.hash);
                return tag;
            } else {
                return {};
            }
        }
                            


    };

};
}
}

#endif
