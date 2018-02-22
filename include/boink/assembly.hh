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
#include "hasher.hh"

#include <deque>
#include <memory>

namespace boink {


typedef std::deque<char> Path;

template <class BaseShifter,
          class GraphType>
class AssemblerMixin : public BaseShifter {

protected:

    std::set<hash_t> seen;

public:
    typedef std::shared_ptr<GraphType> GraphPtr;
    GraphPtr graph;

    AssemblerMixin(GraphPtr graph, BaseShifter const& shifter) :
        BaseShifter(shifter),
        graph(graph) {

    }

    void clear_seen() {
        seen.clear();
    }

    uint8_t degree_left() {
        auto neighbors = this->shift_left();
        return neighbors.size();
    }

    uint8_t degree_right() {
        auto neighbors = this->shift_right();
        return neighbors.size();
    }

    uint8_t degree() {
        return degree_left() + degree_right();
    }

    bool get_left(shift_t& result) {
        std::vector<shift_t> neighbors = this->shift_left();
        if (check_neighbors(neighbors, result)) {
            this->shift_left(result.symbol);
            return true;
        } else {
            return false;
        }
    }

    bool get_right(shift_t& result) {
        std::vector<shift_t> neighbors = this->shift_right();
        if (check_neighbors(neighbors, result)) {
            this->shift_right(result.symbol);
            return true;
        } else {
            return false;
        }
    }

    bool check_neighbors(vector<shift_t> neighbors, shift_t result);

    void assemble_left(const std::string& seed,
                       Path& path) {
        this->reset(seed);
        assemble_left(path);
    }

    void assemble_left(Path& path);

    void assemble_right(const std::string& seed,
                        Path& path) {
        this->reset(seed);
        assemble_right(path);
    }

    void assemble_right(Path& path);

    void assemble(const std::string& seed,
                  Path& path) {
        this->reset(seed);
        assemble(path);
    }
                  
    void assemble(Path& path) {
        assemble_left(path);
        this->get_cursor(path);
        assemble_right(path);
    }
};


template<typename BaseShifter, typename GraphType>
AssemblerMixin<BaseShifter, GraphType> make_assembler(BaseShifter const& shifter) {
    return AssemblerMixin<BaseShifter, GraphType>(shifter);
}

};

#endif
