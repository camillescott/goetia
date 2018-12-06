/* cdbg/metrics.hh -- tracking counters for cDBG stuff
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_CDBG_METRICS_HH
#define BOINK_CDBG_METRICS_HH

#include <string>
#include <iostream>
#include <sstream>

#include <prometheus/registry.h>

#include "boink/cdbg/cdbg_types.hh"

namespace boink {
namespace cdbg {


struct cDBGMetrics {

private:

    std::shared_ptr<prometheus::Registry>        pr_registry;
    prometheus::Family<prometheus::Gauge>&   node_gauge_family;
    prometheus::Family<prometheus::Counter>& op_counter_family;

public:

    prometheus::Gauge&                      n_full;
    prometheus::Gauge&                      n_tips;
    prometheus::Gauge&                      n_islands;
    prometheus::Gauge&                      n_trivial;
    prometheus::Gauge&                      n_circular;
    prometheus::Gauge&                      n_loops;
    prometheus::Gauge&                      n_dnodes;
    prometheus::Gauge&                      n_unodes;

    prometheus::Counter&                    n_splits;
    prometheus::Counter&                    n_merges;
    prometheus::Counter&                    n_extends;
    prometheus::Counter&                    n_clips;
    prometheus::Counter&                    n_deletes;
    prometheus::Counter&                    n_circular_merges;
    
    cDBGMetrics(std::shared_ptr<prometheus::Registry> registry) :
        pr_registry       (registry),
        node_gauge_family (prometheus::BuildGauge()
                                       .Name("boink_cdbg_nodes_current_total")
                                       .Register(*pr_registry)),
        n_full            (node_gauge_family.Add({{"node_type", "full_unode"}})),
        n_tips            (node_gauge_family.Add({{"node_type", "tip_unode"}})),
        n_islands         (node_gauge_family.Add({{"node_type", "island_unode"}})),
        n_trivial         (node_gauge_family.Add({{"node_type", "trivial_unode"}})),
        n_circular        (node_gauge_family.Add({{"node_type", "circular_unode"}})),
        n_loops           (node_gauge_family.Add({{"node_type", "loop_unode"}})),
        n_dnodes          (node_gauge_family.Add({{"node_type", "decision_node"}})),
        n_unodes          (node_gauge_family.Add({{"node_type", "unitig_node"}})),

        op_counter_family (prometheus::BuildCounter()
                                       .Name("boink_cdbg_update_ops_total")
                                       .Register(*pr_registry)),
        n_splits          (op_counter_family.Add({{"update_op", "split"}})),
        n_merges          (op_counter_family.Add({{"update_op", "merge"}})),
        n_extends         (op_counter_family.Add({{"update_op", "extend"}})),
        n_clips           (op_counter_family.Add({{"update_op", "clip"}})),
        n_deletes         (op_counter_family.Add({{"update_op", "delete"}})),
        n_circular_merges (op_counter_family.Add({{"update_op", "merge_circular"}}))
    {
    }

    void increment_cdbg_node(node_meta_t meta) {
        get_cdbg_node_gauage(meta).Increment();
    }

    void decrement_cdbg_node(node_meta_t meta) {
        get_cdbg_node_gauage(meta).Decrement();
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
    prometheus::Gauge& get_cdbg_node_gauage(node_meta_t meta) {
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
    o << c.n_full.Value() << "," 
        << c.n_tips.Value() << "," 
        << c.n_islands.Value() << ","
        << c.n_circular.Value() << ","
        << c.n_trivial.Value() << ","
        << c.n_loops.Value()   << ","
        << c.n_dnodes.Value() << ","
        << c.n_unodes.Value() << ","
        << c.n_splits.Value() << ","
        << c.n_merges.Value() << ","
        << c.n_extends.Value() << ","
        << c.n_clips.Value() << ","
        << c.n_deletes.Value();
    return o;
}

}
}

#endif
