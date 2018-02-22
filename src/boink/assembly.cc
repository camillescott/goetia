#include "boink/assembly.hh"

using namespace boink;

namespace boink {

template <class BaseShifter,
          class GraphType>
bool AssemblerMixin<BaseShifter, GraphType>
::check_neighbors(vector<shift_t> neighbors, shift_t result) {
    
    uint8_t n_found = 0;
    for (auto neighbor : neighbors) {
        if(this->graph->get(neighbor.hash)) {
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


template <class BaseShifter,
          class GraphType>
void AssemblerMixin<BaseShifter, GraphType>
::assemble_right(Path& path) {
    
    if (!graph().get(this->hash())) {
        return;
    }
    shift_t next;
    while (get_right(next) && !seen.count(next.hash)) {
        path.push_back(next.symbol);
        seen.insert(next.hash);
    }
}


template <class BaseShifter,
          class GraphType>
void AssemblerMixin<BaseShifter, GraphType>
::assemble_left(Path& path) {

    if (!graph().get(this->hash())) {
        return;
    }
    shift_t next;
    while (get_left(next) && !seen.count(next.hash)) {
        path.push_front(next.symbol);
        seen.insert(next.hash);
    }
}

};
