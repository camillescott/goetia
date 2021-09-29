/**
 * (c) Camille Scott, 2019
 * File   : breadth_first.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.11.2019
 */

#ifndef GOETIA_BFS_HH
#define GOETIA_BFS_HH

#include "goetia/traversal/unitig_walker.hh"

namespace goetia {
namespace traversal {

template <class GraphType>
struct BreadthFirst : public Traverse<GraphType> {

    using graph_type = GraphType;
    using shifter_type = typename GraphType::shifter_type;
    using hash_type = typename shifter_type::hash_type;
    using shift_type = typename shifter_type::shift_type;
    using kmer_type = typename shifter_type::kmer_type;

    class Search : public Traverse<GraphType>::dBG {

        template<typename Arg1,
                 typename... Args>
        explicit Search(Arg1 arg1, Args&&... args)
            : shifter_type(arg1, std::forward<Args>(args)...)
        {
        }

        explicit Search(shifter_type const& shifter)
            : shifter_type(shifter)
        {
        }

        template<typename... Args>
        static std::shared_ptr<Search> build(Args&&... args) {
            return std::make_shared<Search>(std::forward<Args>(args)...);
        }

        static std::shared_ptr<Search> build(shifter_type const& shifter) {
            return std::make_shared<Search>(shifter);
        }

        //size_t operator() (graph_type * graph


    };

};

}
}

#endif
