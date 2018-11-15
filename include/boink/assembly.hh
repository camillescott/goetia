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

    std::set<hash_t> seen;

public:

    using BaseShifter = typename GraphType::shifter_type;
    using BaseShifter::_K;

    typedef GraphType graph_type;
    typedef BaseShifter shifter_type;

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

    uint8_t degree_left(std::set<hash_t>& extras) {
        auto neighbors = this->gather_left();
        return count_nodes(neighbors, extras);
    }

    uint8_t degree_right(std::set<hash_t>& extras) {
        auto neighbors = this->gather_right();
        return count_nodes(neighbors, extras);
    }

    uint8_t degree(std::set<hash_t>& extras) {
        return degree_left(extras) + degree_right(extras);
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

    uint8_t count_nodes(const vector<shift_t>& nodes,
                        set<hash_t>& extras) {
        uint8_t n_found = 0;
        for (auto node: nodes) {
            if(this->graph->get(node.hash) ||
               extras.count(node.hash)) {
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

    uint8_t reduce_nodes(const vector<shift_t>& nodes,
                         shift_t& result,
                         std::set<hash_t>& extra) {
        uint8_t n_found = 0;
        for (auto node : nodes) {
            //pdebug("check " << neighbor.hash << " " << neighbor.symbol);
            if(this->graph->get(node.hash) ||
               extra.count(node.hash)) {
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

    vector<shift_t> filter_nodes(const vector<shift_t>& nodes,
                                 std::set<hash_t>& extra) {
        vector<shift_t> result;
        for (auto node : nodes) {
            if (this->graph->get(node.hash) ||
                extra.count(node.hash)) {
                result.push_back(node);
            }
        }
        return result;
    }

    vector<kmer_t> find_left_kmers() {
        auto root = this->get_cursor();
        auto filtered = filter_nodes(this->gather_left());
        return graph->build_left_kmers(filtered, root);
    }

    vector<kmer_t> find_right_kmers() {
        auto root = this->get_cursor();
        auto filtered = filter_nodes(this->gather_right());
        return graph->build_right_kmers(filtered, root);
    }

    vector<kmer_t> find_left_kmers(std::set<hash_t>& extras) {
        auto root = this->get_cursor();
        auto filtered = filter_nodes(this->gather_left(), extras);
        return graph->build_left_kmers(filtered, root);
    }

    vector<kmer_t> find_right_kmers(std::set<hash_t>& extras) {
        auto root = this->get_cursor();
        auto filtered = filter_nodes(this->gather_right(), extras);
        return graph->build_right_kmers(filtered, root);
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


template<class GraphType>
class CompactorMixin : public AssemblerMixin<GraphType> {

public:

    using BaseShifter = typename GraphType::shifter_type;
    using BaseShifter::_K;
    using AssemblerType = AssemblerMixin<GraphType>;
    using AssemblerType::seen;
    using AssemblerType::get_left;
    using AssemblerType::get_right;
    using AssemblerType::degree_left;
    using AssemblerType::degree_right;
    using AssemblerType::count_nodes;
    using AssemblerType::filter_nodes;
    using AssemblerType::find_left_kmers;
    using AssemblerType::find_right_kmers;
    using AssemblerType::gather_left;
    using AssemblerType::gather_right;

    typedef AssemblerType assembler_type;
    typedef BaseShifter shifter_type;

    CompactorMixin(GraphType * graph, BaseShifter const& shifter) :
        AssemblerType(graph, shifter)
    {
    }

    CompactorMixin(GraphType * graph) :
        AssemblerType(graph)
    {
    }

    string compactify(const string& seed) {
        Path path;
        this->set_cursor(seed);
        this->get_cursor(path);
        hash_t end_hash;
        compactify_left(path, end_hash);
        this->set_cursor(seed);
        compactify_right(path, end_hash);

        return this->to_string(path);
    }

    void compactify_right(Path& path, hash_t& end_hash, set<hash_t>& mask) {
        end_hash = this->get();
        this->seen.clear();
        this->seen.insert(this->get());
        
        shift_t next;
        uint8_t n_right;
        while (1) {
            if (degree_left() > 1) {
                path.pop_back();
                return;
            }

            n_right = this->reduce_nodes(this->gather_right(),
                                         next);
            if (n_right > 1) {
                path.pop_back();
                return;
            }

            if (n_right == 0) {
                end_hash = this->get();
                return;
            }

            if (this->seen.count(next.hash) ||
                mask.count(next.hash)) {
                end_hash = this->get();
                return;
            }

            end_hash = this->get();
            this->shift_right(next.symbol);
            path.push_back(next.symbol);
            this->seen.insert(next.hash);
        }
    }

    void compactify_left(Path& path, hash_t& end_hash, set<hash_t>& mask) {
        end_hash = this->get();
        this->seen.clear();
        this->seen.insert(this->get());

        shift_t next;
        uint8_t n_left;
        while (1) {
            if (degree_right() > 1) {
                pdebug("Stop: reverse d-node");
                path.pop_front();
                return;
            }

            n_left = this->reduce_nodes(this->gather_left(),
                                        next);
            if (n_left > 1) {
                pdebug("Stop: forward d-node");
                path.pop_front();
                return;
            }

            if (n_left == 0) {
                pdebug("Stop: no neighbors.");
                end_hash = this->get();
                return;
            }

            if (this->seen.count(next.hash) ||
                mask.count(next.hash)) {
                pdebug("Stop: hit seen k-mer or masked k-mer.");
                end_hash = this->get();
                return;
            }
            
            end_hash = this->get();
            this->shift_left(next.symbol);
            path.push_front(next.symbol);
            this->seen.insert(next.hash);
        }
    }

    bool is_decision_kmer(const string& node,
                          uint8_t& degree) {
        this->set_cursor(node);
        return is_decision_kmer(degree);
    }

    bool is_decision_kmer(const string& node) {
        this->set_cursor(node);
        return this->degree_left() > 1 || this->degree_right() > 1;
    }

    bool is_decision_kmer(uint8_t& degree) {
        uint8_t ldegree, rdegree;
        ldegree = this->degree_left();
        rdegree = this->degree_right();
        degree = ldegree + rdegree;
        return ldegree > 1 || rdegree > 1;
    }

    void find_decision_kmers(const string& sequence,
                             vector<uint32_t>& decision_positions,
                             HashVector& decision_hashes,
                             vector<NeighborBundle>& decision_neighbors) {

        KmerIterator<AssemblerType> iter(sequence, this);
        size_t pos = 0;
        while(!iter.done()) {
            hash_t h = iter.next();
            NeighborBundle neighbors;
            if (get_decision_neighbors(iter.shifter,
                                       neighbors)) {

                decision_neighbors.push_back(neighbors);
                decision_positions.push_back(pos);
                decision_hashes.push_back(h);
            }
        
            ++pos;
       }
    }

    bool get_decision_neighbors(const string& root,
                                NeighborBundle& result,
                                set<hash_t>& union_nodes) {
        return get_decision_neighbors(this, root, result, union_nodes);
    }

    bool get_decision_neighbors(NeighborBundle& result, set<hash_t>& union_nodes) {
        return get_decision_neighbors(this, result, union_nodes);
    }

    bool get_decision_neighbors(AssemblerType* shifter,
                                const string& root,
                                NeighborBundle& result,
                                set<hash_t>& union_nodes) {
        shifter->set_cursor(root);
        return get_decision_neighbors(shifter, result, union_nodes);
    }

    bool get_decision_neighbors(AssemblerType* shifter,
                                NeighborBundle& result,
                                set<hash_t>& union_nodes) {

        auto left_kmers = shifter->find_left_kmers(union_nodes);
        auto right_kmers = shifter->find_right_kmers(union_nodes);
        
        if (left_kmers.size() > 1 || right_kmers.size() > 1) {
            result = std::make_pair(left_kmers, right_kmers);
            return true;
        } else {
            return false;
        }
    }

    bool get_decision_neighbors(const string& root,
                                NeighborBundle& result) {
        return get_decision_neighbors(this, root, result);
    }

    bool get_decision_neighbors(NeighborBundle& result) {
        return get_decision_neighbors(this, result);
    }

    bool get_decision_neighbors(AssemblerType* shifter,
                                const string& root,
                                NeighborBundle& result) {
        shifter->set_cursor(root);
        return get_decision_neighbors(shifter, result);
    }

    bool get_decision_neighbors(AssemblerType* shifter,
                                NeighborBundle& result) {

        auto left_kmers = shifter->find_left_kmers();
        auto right_kmers = shifter->find_right_kmers();
        
        if (left_kmers.size() > 1 || right_kmers.size() > 1) {
            result = std::make_pair(left_kmers, right_kmers);
            return true;
        } else {
            return false;
        }
    }
};

}

#undef pdebug
#endif
