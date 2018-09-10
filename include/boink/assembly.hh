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
typedef std::pair<kmer_t, NeighborBundle> DecisionKmer;

struct DecisionKmerHash {
    // Just uses the hash of the kmer_t from the pair
    inline std::size_t operator()(const DecisionKmer& n) const {
        return n.first.hash;
    }
};

struct  DecisionKmerComp {
    bool operator()(const DecisionKmer& a, const DecisionKmer& b) {
        return a.first.hash < b.first.hash;
    }
};

typedef std::set<DecisionKmer, DecisionKmerComp> DecisionKmerSet;
typedef std::unordered_set<DecisionKmer, DecisionKmerHash> DecisionKmerUSet;


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

    uint8_t get_left(shift_t& result) {
        std::vector<shift_t> neighbors = this->gather_left();
        auto n_left = reduce_nodes(neighbors, result);

        if (n_left == 1) {
            this->shift_left(result.symbol);
        }
        return n_left;
    }

    uint8_t get_right(shift_t& result) {
        std::vector<shift_t> neighbors = this->gather_right();
        auto n_right = reduce_nodes(neighbors, result);

        if (n_right == 1) {
            this->shift_right(result.symbol);
        }
        return n_right;
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

    uint8_t reduce_nodes(const vector<shift_t>& nodes,
                      shift_t& result) {
        uint8_t n_found = 0;
        for (auto node : nodes) {
            //pdebug("check " << neighbor.hash << " " << neighbor.symbol);
            if(this->graph->get(node.hash)) {
                //pdebug("found " << neighbor.hash);
                ++n_found;
                if (n_found > 1) {
                    return n_found;
                }
                result = node;
            }
        }
        return n_found;
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

    vector<kmer_t> build_left_kmers(const vector<shift_t>& nodes,
                                    const string& root) {
        KmerVector kmers;
        for (auto neighbor : nodes) {
            kmers.push_back(kmer_t(neighbor.hash,
                                   neighbor.symbol
                                   + root.substr(0, this->_K-1)));
        }
        return kmers;
    }

    vector<kmer_t> build_right_kmers(const vector<shift_t>& nodes,
                                     const string& root) {
        KmerVector kmers;
        for (auto neighbor : nodes) {
            kmers.push_back(kmer_t(neighbor.hash,
                                   root.substr(1, this->_K-1)
                                   + neighbor.symbol));
        }

        return kmers;
    }

    vector<kmer_t> find_left_kmers() {
        auto root = this->get_cursor();
        auto filtered = filter_nodes(this->gather_left());
        return build_left_kmers(filtered, root);
    }

    vector<kmer_t> find_right_kmers() {
        auto root = this->get_cursor();
        auto filtered = filter_nodes(this->gather_right());
        return build_right_kmers(filtered, root);
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
