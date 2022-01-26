/* cdbg/metrics.hh -- tracking counters for cDBG stuff
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef GOETIA_CDBG_METRICS_HH
#define GOETIA_CDBG_METRICS_HH

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include "goetia/metrics.hh"
#include "goetia/cdbg/cdbg_types.hh"

namespace goetia {


class cDBGMetrics {

public:

    Gauge                    n_full;
    Gauge                    n_tips;
    Gauge                    n_islands;
    Gauge                    n_trivial;
    Gauge                    n_circular;
    Gauge                    n_loops;
    Gauge                    n_dnodes;
    Gauge                    n_unodes;


    Gauge                    n_splits;
    Gauge                    n_merges;
    Gauge                    n_extends;
    Gauge                    n_clips;
    Gauge                    n_deletes;
    Gauge                    n_circular_merges;
    
    cDBGMetrics() :
        n_full            {"node_type", "full_unode"},
        n_tips            {"node_type", "tip_unode"},
        n_islands         {"node_type", "island_unode"},
        n_trivial         {"node_type", "trivial_unode"},
        n_circular        {"node_type", "circular_unode"},
        n_loops           {"node_type", "loop_unode"},
        n_dnodes          {"node_type", "decision_node"},
        n_unodes          {"node_type", "unitig_node"},

        n_splits          {"update_op", "split"},
        n_merges          {"update_op", "merge"},
        n_extends         {"update_op", "extend"},
        n_clips           {"update_op", "clip"},
        n_deletes         {"update_op", "delete"},
        n_circular_merges {"update_op", "merge_circular"}
    {
    }

    void increment_cdbg_node(node_meta_t meta) {
        get_cdbg_node_gauage(meta)++;
    }

    void decrement_cdbg_node(node_meta_t meta) {
        get_cdbg_node_gauage(meta)--;
    }

    std::string repr() const {
        std::ostringstream os;
        os << *this;
        return os.str();
    }

    std::string header() const {
        return std::string("n_full,n_tips,n_islands,n_circular,n_trivial,n_loops,"
                           "n_dnodes,n_unodes,n_splits,n_merges,n_extends,n_clips,n_deletes");
    }

    friend std::ostream& operator<<(std::ostream& o, const cDBGMetrics& c);

private:
    Gauge& get_cdbg_node_gauage(node_meta_t meta) {
        switch(meta) {
            case FULL:
                return n_full;
            case TIP:
                return n_tips;
            case ISLAND:
                return n_islands;
            case CIRCULAR:
                return n_circular;
            case TRIVIAL:
                return n_trivial;
            case LOOP:
                return n_loops;
            case DECISION:
                return n_dnodes;
        }
    }
};


inline std::ostream& operator<<(std::ostream& o, const cDBGMetrics& c) {
    o << c.n_full << "," 
      << c.n_tips << "," 
      << c.n_islands << ","
      << c.n_circular << ","
      << c.n_trivial << ","
      << c.n_loops << ","
      << c.n_dnodes << ","
      << c.n_unodes << ","
      << c.n_splits << ","
      << c.n_merges << ","
      << c.n_extends << ","
      << c.n_clips << ","
      << c.n_deletes;
    return o;
}

}

#endif
