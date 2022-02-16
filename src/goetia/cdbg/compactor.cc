/**
 * (c) Camille Scott, 2019
 * File   : compactor.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.01.2020
 */

#include "goetia/cdbg/compactor.hh"
#include "goetia/storage/storage_types.hh"
#include "goetia/hashing/rollinghashshifter.hh"

namespace goetia {

# ifdef DEBUG_CPTR
#   define pdebug(x) do { std::cerr << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif

template <template <class, class> class GraphType,
            class StorageType,
            class ShifterType>
void
StreamingCompactor<GraphType<StorageType, ShifterType>>::
Compactor::_induce_decision_nodes(std::deque<DecisionKmer>& induced_decision_kmers,
                                  std::set<hash_type>& new_kmers) {

    std::set<hash_type> induced_decision_kmer_hashes;

    pdebug("Perform induction on " << induced_decision_kmers.size() <<
           " new decision k-mers");
    bool found = false;
    for (auto dkmer : induced_decision_kmers) {
        _build_dnode(dkmer.first);
        induced_decision_kmer_hashes.insert(dkmer.first.value());
    }

    std::set<hash_type> processed;
    size_t n_attempts = 0;
    size_t max_attempts = 4 * induced_decision_kmer_hashes.size();
    while (induced_decision_kmers.size() > 0) {
        pdebug(induced_decision_kmers.size() << " more splits to attempt...");
        n_attempts++;

        DecisionKmer d_kmer(induced_decision_kmers.front());
        induced_decision_kmers.pop_front();

        if (processed.count(d_kmer.first.value())) {
            pdebug("Processed " << d_kmer.first << " already");
            continue;
        }

        if (_try_split_unode(d_kmer.first,
                             d_kmer.second, 
                             new_kmers,
                             induced_decision_kmer_hashes,
                             processed)) {
            pdebug("Split successful on " << d_kmer.first);
            processed.insert(d_kmer.first.value());
        } else {
            induced_decision_kmers.push_back(d_kmer);
        }
        
        if (n_attempts > max_attempts) {
            throw GoetiaException("Stuck in split attempt loop, failing.");
        }
    }
}



template <template <class, class> class GraphType,
            class StorageType,
            class ShifterType>
bool
StreamingCompactor<GraphType<StorageType, ShifterType>>::
Compactor::_try_split_unode(kmer_type root,
                            neighbor_pair_type& neighbors,
                            std::set<hash_type>& new_kmers,
                            std::set<hash_type>& induced_decision_kmer_hashes,
                            std::set<hash_type>& processed) {
    pdebug("Attempt unitig split from " << root);

    UnitigNode * unode_to_split;
    
    if ((unode_to_split = cdbg->query_unode_end(root)) != nullptr) {
        // special case: induced an end k-mer, just have to trim the u-node,
        // no need to create a new one

        if (unode_to_split->meta() == TRIVIAL ||
            unode_to_split->sequence.size() == this->K) {
            pdebug("Induced a trivial u-node, delete it.");
            cdbg->delete_unode(unode_to_split);
            return true;
        }

        if (unode_to_split->meta() == CIRCULAR) {
            auto unitig = unode_to_split->sequence;
            cdbg->split_unode(unode_to_split->node_id,
                              0,
                              root.kmer,
                              dbg->hash(unitig.substr(unitig.size() - this->K)),
                              dbg->hash(unitig.c_str() + 1));
            return true;
        }

        hash_type new_end;
        bool clip_from;
        if (root.value() == unode_to_split->left_end().value()) {
            new_end = dbg->hash(unode_to_split->sequence.c_str() + 1);
            clip_from = DIR_LEFT;
        } else {
            new_end = dbg->hash(unode_to_split->sequence.c_str()
                                 + unode_to_split->sequence.size()
                                 - this->K - 1);
            clip_from = DIR_RIGHT;
        }

        cdbg->clip_unode(clip_from,
                         root,
                         new_end);

        return true;
    }

    std::vector<kmer_type> lfiltered;
    std::copy_if(neighbors.first.begin(),
                 neighbors.first.end(),
                 std::back_inserter(lfiltered),
                 [&] (kmer_type neighbor) { return
                    !new_kmers.count(neighbor) &&
                    !processed.count(neighbor);
                 });

    std::vector<kmer_type> rfiltered;
    std::copy_if(neighbors.second.begin(),
                 neighbors.second.end(),
                 std::back_inserter(rfiltered),
                 [&] (kmer_type neighbor) { return
                    !new_kmers.count(neighbor) &&
                    !processed.count(neighbor);
                 });
    // EDGE CASE: split k-mer has no unitig to split because of previous processed
    // nodes or because neighbor is a properly oriented unitig end

    pdebug(lfiltered.size() << " left, " << rfiltered.size() << " right");

    //auto masked = Masked<StorageType, ShifterType, std::set<hash_type>>(*dbg, processed);
    WalkStopper<hash_type> stopper_func(processed);

    if (lfiltered.size()) {
        // size should always be 1 here
        pdebug("Found a valid left neighbor, search this way... ("
               << lfiltered.size() << " in filtered set, should always be 1.)");
        auto start = lfiltered.back();
        //auto walk = masked.walk_left(start.kmer);
        walk_type<DIR_LEFT> walk = dbg->walk_left(start.kmer, stopper_func);

        hash_type end_hash = walk.tail();

        if (walk.end_state != State::STOP_SEEN) {
            dbg->clear_seen();

            if (walk.end_state == State::DECISION_FWD) {
                //auto step = masked.step_right();
                auto step = dbg->step_right();
                pdebug("stopped on DECISION_FWD " << step.first);
                if (step.first == State::STEP) {
                    pdebug("pop off path");
                    end_hash = step.second.front();
                    walk.path.pop_back();
                }
            }
            unode_to_split = cdbg->query_unode_end(end_hash);

            size_t split_point = walk.path.size() + 1;
            hash_type left_unode_new_right = start;

            pdebug("split point is " << split_point <<
                    " new_right is " << left_unode_new_right
                   << " root was " << root);
            assert(unode_to_split != nullptr);

            hash_type right_unode_new_left = dbg->hash(unode_to_split->sequence.c_str() + 
                                                        split_point + 1,
                                                        this->K);

            if (rfiltered.size()) {
                if (right_unode_new_left.value() != rfiltered.back().value()) {
                  std::cout << "unode to split: " << *unode_to_split << std::endl;
                  std::cout << "walk end hash: " << end_hash << std::endl;
                  std::cout << "walk state: " << walk.end_state << std::endl;
                  std::cout << "right unode's new left end: " << right_unode_new_left.value() << std::endl;
                  std::cout << "right unode's new left kmer: " << unode_to_split->sequence.substr(split_point + 1, this->K) << std::endl;
                  std::cout << "rfiltered size: " << rfiltered.size() << std::endl;
                  std::cout << "start: " << start << std::endl;
                  std::cout << "root: " << root << std::endl;
                  std::cout << "start in induced: " << induced_decision_kmer_hashes.count(start) << std::endl;
                }
                assert(right_unode_new_left.value() == rfiltered.back().value());
            }

            cdbg->split_unode(unode_to_split->node_id,
                              split_point,
                              root.kmer,
                              left_unode_new_right,
                              right_unode_new_left);

            return true;
        } else {
            pdebug("Unitig to split is a loop, traversed " << dbg->seen.size());
            for (auto hash : dbg->seen) {
                unode_to_split = cdbg->query_unode_end(hash);
                if (unode_to_split != nullptr) break;
            }
            if (unode_to_split != nullptr) {
                pdebug("u-node to split: " << *unode_to_split);
                cdbg->split_unode(unode_to_split->node_id,
                                  0,
                                  root.kmer,
                                  lfiltered.back(),
                                  rfiltered.back());
                return true;
            }
        }
    }

    if (rfiltered.size()) {
        // size should always be 1 here
        pdebug("Found a valid right neighbor, search this way... ("
               << rfiltered.size() << " in filtered set, should be 1.");
        auto start = rfiltered.back();
        //auto walk = masked.walk_right(start.kmer);
        walk_type<DIR_RIGHT> walk = dbg->walk_right(start.kmer, stopper_func);
        hash_type end_hash = walk.tail();

        if (walk.end_state != State::STOP_SEEN) {
            dbg->clear_seen();

            if (walk.end_state == State::DECISION_FWD) {

                //auto step = masked.step_left();
                auto step = dbg->step_left();
                if (step.first == State::STEP) {
                    end_hash = step.second.front();
                    walk.path.pop_back();
                }
            }

            unode_to_split = cdbg->query_unode_end(end_hash);
            size_t split_point = unode_to_split->sequence.size()
                                                 - walk.path.size()
                                                 - this->K 
                                                 - 1;
            hash_type new_right = dbg->hash(unode_to_split->sequence.c_str() + 
                                             split_point - 1,
                                             this->K);
            hash_type new_left = start;
            if (lfiltered.size()) {
                assert(lfiltered.back().value() == new_right.value());
            }

            cdbg->split_unode(unode_to_split->node_id,
                              split_point,
                              root.kmer,
                              new_right,
                              new_left);

            return true;
        } else {
            pdebug("Unitig to split is a loop, traversed " << dbg->seen.size());
            for (auto hash : dbg->seen) {
                unode_to_split = cdbg->query_unode_end(hash);
                if (unode_to_split != nullptr) break;
            }
            if (unode_to_split != nullptr) {
                pdebug("u-node to split: " << *unode_to_split);
                cdbg->split_unode(unode_to_split->node_id,
                                  0,
                                  root.kmer,
                                  lfiltered.back(),
                                  rfiltered.back());
                return true;
            }
        }
    }

    // failed to find a split point. must be flanked by induced d-nodes on either
    // side before a tag or unode end: aka same unode split by three or more
    // induced d-nodes
    pdebug("Attempt to split unode failed from " << root);

    return false;
}


}


template class goetia::StreamingCompactor<goetia::dBG<goetia::SparseppSetStorage, goetia::FwdLemireShifter>>;
template class goetia::StreamingCompactor<goetia::dBG<goetia::PHMapStorage, goetia::FwdLemireShifter>>;
// template class goetia::StreamingCompactor<goetia::dBG<goetia::BitStorage, goetia::FwdLemireShifter>>;
// template class goetia::StreamingCompactor<goetia::dBG<goetia::ByteStorage, goetia::FwdLemireShifter>>;
// template class goetia::StreamingCompactor<goetia::dBG<goetia::NibbleStorage, goetia::FwdLemireShifter>>;
// template class goetia::StreamingCompactor<goetia::dBG<goetia::QFStorage, goetia::FwdLemireShifter>>;

