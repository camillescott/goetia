/* assembly.hh -- boink traversal and assembly
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef ASSEMBLY_HH
#define ASSEMBLY_HH

#include "hashing.hh"

#include <deque>


# ifdef DEBUG_ASMLY
#   define pdebug(x) do { std::cerr << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif


namespace boink {


typedef std::deque<char> Path;
typedef std::vector<std::string> StringVector;
typedef std::vector<kmer_t> KmerVector;
typedef std::pair<KmerVector, KmerVector> NeighborBundle;


template<class GraphType>
class AssemblerMixin : public GraphType::shifter_type {

protected:

    using BaseShifter = typename GraphType::shifter_type;
    using BaseShifter::_K;
    std::set<hash_t> seen;

public:

    GraphType * graph;

    AssemblerMixin(GraphType * graph, BaseShifter const& shifter) :
        BaseShifter(shifter),
        graph(graph) {

    }

    AssemblerMixin(GraphType * graph) :
        BaseShifter(graph->K()),
        graph(graph) {
    }

    void clear_seen() {
        seen.clear();
    }

    uint8_t degree_left() {
        auto neighbors = this->gather_left();
        return count_nodes(neighbors);
    }

    uint8_t degree_right() {
        auto neighbors = this->gather_right();
        return count_nodes(neighbors);
    }

    uint8_t degree() {
        return degree_left() + degree_right();
    }

    bool get_left(shift_t& result) {
        std::vector<shift_t> neighbors = this->gather_left();
        if (reduce_nodes(neighbors, result)) {
            this->shift_left(result.symbol);
            return true;
        } else {
            return false;
        }
    }

    bool get_right(shift_t& result) {
        std::vector<shift_t> neighbors = this->gather_right();
        if (reduce_nodes(neighbors, result)) {
            this->shift_right(result.symbol);
            return true;
        } else {
            return false;
        }
    }

    uint8_t count_nodes(const vector<shift_t>& nodes) {
        uint8_t n_found = 0;
        for (auto node: nodes) {
            if(this->graph->get(node.hash)) {
                ++n_found;
            }
        }
        return n_found;
    }

    bool reduce_nodes(const vector<shift_t>& nodes,
                      shift_t& result) {
        uint8_t n_found = 0;
        for (auto node : nodes) {
            //pdebug("check " << neighbor.hash << " " << neighbor.symbol);
            if(this->graph->get(node.hash)) {
                //pdebug("found " << neighbor.hash);
                ++n_found;
                if (n_found > 1) {
                    return false;
                }
                result = node;
            }
        }
        if (n_found == 0) {
            return false;
        } else {
            return true;
        }
    }

    vector<shift_t> filter_nodes(const vector<shift_t>& nodes) {
        vector<shift_t> result;
        for (auto node : nodes) {
            if (this->graph->get(node.hash)) {
                result.push_back(node);
            }
        }
        return result;
    }

    void assemble_left(const std::string& seed,
                       Path& path) {

        this->set_cursor(seed);
        if (!this->graph->get(this->get())) {
            return;
        } 
        this->get_cursor(path);
        assemble_left(path);
    }

    void assemble_left(Path& path) {
        seen.insert(this->get());

        shift_t next;
        while (get_left(next) && !seen.count(next.hash)) {
            path.push_front(next.symbol);
            seen.insert(next.hash);
        }
    }

    void assemble_right(const std::string& seed,
                        Path& path) {

        this->set_cursor(seed);
        if (!this->graph->get(this->get())) {
            return;
        }
        this->get_cursor(path);
        assemble_right(path);
    }

    void assemble_right(Path& path) {
        seen.insert(this->get());

        shift_t next;
        while (get_right(next) && !seen.count(next.hash)) {
            path.push_back(next.symbol);
            seen.insert(next.hash);
        }
    }

    void assemble(const std::string& seed,
                  Path& path) {

        this->set_cursor(seed);
        if (!this->graph->get(this->hash(seed))) {
            return;
        }
        this->get_cursor(path);
        assemble_left(path);
        this->set_cursor(seed);
        assemble_right(path);
    }
                  
    string to_string(Path& path) {
        return string(path.begin(), path.end());
    }
};


template<class dBGType>
AssemblerMixin<dBGType> make_assembler(dBGType * graph) {
    return AssemblerMixin<dBGType>(graph);
}

}

#undef pdebug
#endif
