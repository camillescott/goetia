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

namespace boink {


typedef std::deque<char> Path;

template <class BaseShifter,
          class GraphType>
class AssemblerMixin {

private:

    BaseShifter& shifter() {
        return *static_cast<BaseShifter*>(this);
    }

public:
    typedef std::shared_ptr<GraphType> GraphPtr;
    GraphPtr graph;

    AssemblerMixin(GraphPtr graph) :
        graph(graph) {

    }


    uint8_t degree_left() {
        std::vector<hash_t> neighbors = shifter().shift_left();
        return neighbors.size();
    }

    uint8_t degree_right() {
        std::vector<hash_t> neighbors = shifter().shift_right();
        return neighbors.size();
    }

    uint8_t degree() {
        return degree_left() + degree_right();
    }

    std::string assemble_left(const std::string& seed,
                              ) {
        reset(seed);
        return assemble_left();

    std::string assemble_left() {
        if (!graph().get(shifter.hash())) {
            return "";
        }
        std::vector<hash_t> neighbors = shifter.shift_left();
        while(neighbors.size() == 1) {
            path_deque.
        }
    }



};

};

#endif
