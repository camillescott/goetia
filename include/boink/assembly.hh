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

#include "dbg.hh"
#include "hashing.hh"

#include <deque>
#include <memory>


# ifdef DEBUG_ASSEMBLY
#   define pdebug(x) do { std::cerr << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif


namespace boink {


typedef std::deque<char> Path;


template<class GraphType>
class AssemblerMixin : public GraphType::shifter_type {

protected:

    std::set<hash_t> seen;

public:

    using BaseShifter = typename GraphType::shifter_type;
    typedef std::shared_ptr<GraphType> GraphPtr;
    GraphPtr graph;

    AssemblerMixin(GraphPtr graph, BaseShifter const& shifter) :
        BaseShifter(shifter),
        graph(graph) {

    }

    AssemblerMixin(GraphPtr graph) :
        BaseShifter(graph->K()),
        graph(graph) {
    }

    void clear_seen() {
        seen.clear();
    }

    uint8_t degree_left() {
        auto neighbors = this->gather_left();
        return neighbors.size();
    }

    uint8_t degree_right() {
        auto neighbors = this->gather_right();
        return neighbors.size();
    }

    uint8_t degree() {
        return degree_left() + degree_right();
    }

    bool get_left(shift_t& result) {
        std::vector<shift_t> neighbors = this->gather_left();
        if (check_neighbors(neighbors, result)) {
            this->shift_left(result.symbol);
            return true;
        } else {
            return false;
        }
    }

    bool get_right(shift_t& result) {
        std::vector<shift_t> neighbors = this->gather_right();
        if (check_neighbors(neighbors, result)) {
            this->shift_right(result.symbol);
            return true;
        } else {
            return false;
        }
    }

    bool check_neighbors(vector<shift_t> neighbors,
                         shift_t& result) {
        uint8_t n_found = 0;
        for (auto neighbor : neighbors) {
            //pdebug("check " << neighbor.hash << " " << neighbor.symbol);
            if(this->graph->get(neighbor.hash)) {
                //pdebug("found " << neighbor.hash);
                ++n_found;
                if (n_found > 1) {
                    return false;
                }
                result = neighbor;
            }
        }
        if (n_found == 0) {
            return false;
        } else {
            return true;
        }
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
AssemblerMixin<dBGType> make_assembler(shared_ptr<dBGType> graph) {
    return AssemblerMixin<dBGType>(graph);
}

}

#endif
