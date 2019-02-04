# boink/tests/test_cdbg.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import itertools
import sys

import pytest
from boink.tests.utils import *

from boink.compactor import display_segment_list, make_streaming_compactor
from boink.prometheus import Instrumentation


@pytest.fixture
def compactor(ksize, graph, graph_type):
    _graph_type, AdapterType = graph_type
    instrumention = Instrumentation('', expose=False)
    compactor = make_streaming_compactor(graph, inst=instrumention)
    return compactor


class TestFindNewSegments:

    @using_ksize(15)
    @using_length(100)
    @pytest.mark.benchmark(group='cdbg-segments')
    def test_fork_core_first(self, ksize, length, graph, compactor, right_fork,
                                   check_fp, benchmark):
        (core, branch), pivot = right_fork()
        check_fp()

        print('INPUTS', core, core[:pivot+1], ' ' * (pivot + 1) + branch, sep='\n\n')

        core_segments = compactor.find_new_segments(core)
        graph.add_sequence(core)
        branch_segments = benchmark(compactor.find_new_segments, 
                                    core[:pivot+1] + branch)
        graph.add_sequence(core[:pivot+1] + branch)

        # should be NULL - SEG - NULL
        print('Core Segments')
        display_segment_list(core_segments)
        for s in core_segments:
            print('segment:', s.sequence)
        # again NULL - SEG - NULL
        print('Branch Segments')
        display_segment_list(branch_segments)
        for s in branch_segments:
            print('segment:', s.sequence)

        assert branch_segments[0].is_null
        assert len(core_segments) == 3
        assert len(branch_segments) == 3
        assert branch_segments[1].sequence == branch
        assert branch_segments[1].start == pivot + 1
        assert branch_segments[1].is_decision_kmer == False
        assert branch_segments[-1].is_null

        assert core_segments[0].is_null
        assert core_segments[1].sequence == core
        assert len(core_segments[1].sequence) == length
        assert core_segments[1].start == 0
        assert core_segments[1].length == length
        assert core_segments[1].is_decision_kmer == False
        assert core_segments[-1].is_null

        for segment in branch_segments:
            if not segment.is_null:
                for kmer in kmers(segment.sequence, ksize):
                    assert graph.left_degree(kmer) < 2
                    assert graph.right_degree(kmer) < 2

    @using_ksize(15)
    @using_length(100)
    @pytest.mark.benchmark(group='cdbg-segments')
    def test_right_decision_split(self, ksize, graph, compactor, right_fork,
                                  check_fp, benchmark):
        (core, branch), pivot = right_fork()
        check_fp()

        print('INPUTS', core, core[:pivot+1], ' ' * (pivot + 1) + branch, sep='\n\n')

        branch_segments = compactor.find_new_segments(branch)
        graph.add_sequence(branch)
        core_segments = benchmark(compactor.find_new_segments, core)
        graph.add_sequence(core)

        print('Branch Segments')
        display_segment_list(branch_segments)
        for s in branch_segments:
            print('segment:' + s.sequence)

        print('Core Segments')
        display_segment_list(core_segments)
        for s in core_segments:
            print('segment:' + s.sequence)

        assert len(branch_segments) == 3
        assert branch_segments[0].is_null
        assert branch_segments[-1].is_null
        branch_segment = branch_segments[1]
        assert branch_segment.sequence == branch
        assert branch_segment.length == len(branch)

        assert len(core_segments) == 5
        assert core_segments[0].is_null
        assert core_segments[1].is_decision_kmer == False
        assert core_segments[2].is_decision_kmer == True
        assert core_segments[3].is_decision_kmer == False
        assert core_segments[-1].is_null

        assert graph.left_degree(core_segments[2].sequence) == 1
        assert graph.right_degree(core_segments[2].sequence) == 2
        assert graph.left_degree(core[pivot:pivot+ksize]) == 1
        assert graph.right_degree(core[pivot:pivot+ksize]) == 2

        assert core[:pivot+ksize-1] == core_segments[1].sequence
        assert core[pivot+1:] == core_segments[3].sequence


    @using_ksize(15)
    @using_length(100)
    @pytest.mark.benchmark(group='cdbg-segments')
    def test_merge_no_decisions(self, ksize, length, graph, compactor,
                                      linear_path, check_fp, benchmark):
        sequence = linear_path()
        left = sequence[:length//2]
        right = sequence[length//2:]
        print(left)
        print(right)
        check_fp()

        left_segments = compactor.find_new_segments(left)
        graph.add_sequence(left)
        right_segments = compactor.find_new_segments(right)
        graph.add_sequence(right)

        assert len(left_segments) == 3
        assert left_segments[0].is_null
        assert left_segments[1].sequence == left
        assert left_segments[1].left_anchor == graph.hash(left[:ksize])
        assert left_segments[1].right_anchor == graph.hash(left[-ksize:])
        assert left_segments[-1].is_null

        assert right_segments[0].is_null
        assert len(right_segments) == 3
        assert right_segments[1].sequence == right
        assert right_segments[-1].is_null

        merged = left + right
        merged_new = left[-ksize+1:] + right[:ksize-1]
        assert sum(graph.get_counts(merged_new)) == 0
        merged_segments = benchmark(compactor.find_new_segments, merged)
        assert len(merged_segments) == 3
        assert merged_segments[0].is_null
        assert merged_segments[-1].is_null
        merged_segment = merged_segments[1]

        assert merged_segment.sequence == merged_new
        assert merged_segment.left_anchor == graph.hash(merged_new[:ksize])
        assert merged_segment.right_anchor == graph.hash(merged_new[-ksize:])


    @using_ksize(15)
    @using_length(100)
    @pytest.mark.benchmark(group='cdbg-segments')
    def test_right_decision_on_end(self, ksize, graph, compactor, right_fork,
                                         check_fp, benchmark):
        (core, branch), pivot = right_fork()
        check_fp()

        # pivot is start position of decision k-mer
        print('INPUTS', core, core[:pivot+1], ' ' * (pivot + 1) + branch, sep='\n\n')
        upper = branch
        lower = core[pivot+1:]
        test = core[:pivot+ksize]

        upper_segments = compactor.find_new_segments(upper)
        graph.add_sequence(upper)
        lower_segments = compactor.find_new_segments(lower)
        graph.add_sequence(lower)

        assert len(upper_segments) == 3
        assert upper_segments[1].sequence == upper
        assert len(lower_segments) == 3
        assert lower_segments[1].sequence == lower

        test_segments = benchmark(compactor.find_new_segments, test)
        display_segment_list(test_segments)

        assert test_segments[0].is_null
        assert len(test_segments) == 4
        assert test_segments[2].is_decision_kmer
        assert len(test_segments[2].sequence) == ksize
        assert len(test_segments[1].sequence) == pivot + ksize - 1
        assert test_segments[1].left_anchor == graph.hash(core[:ksize])
        assert test_segments[1].right_anchor == graph.hash(core[pivot-1:pivot+ksize-1])
        assert test_segments[-1].is_null


    @using_ksize(15)
    @using_length(100)
    @pytest.mark.benchmark(group='cdbg-segments')
    def test_left_decision_on_end(self, ksize, length, graph, compactor, left_fork,
                                        check_fp, benchmark):
        (core, branch), pivot = left_fork()
        check_fp()
        # pivot is start position of decision k-mer
        print('INPUTS', core, core[:pivot+1], ' ' * (pivot + 1) + branch, sep='\n\n')
        upper = branch
        lower = core[:pivot+ksize-1]
        test = core[pivot:]

        upper_segments = compactor.find_new_segments(upper)
        graph.add_sequence(upper)
        lower_segments = compactor.find_new_segments(lower)
        graph.add_sequence(lower)

        assert len(upper_segments) == 3
        assert upper_segments[1].sequence == upper
        assert len(lower_segments) == 3
        assert lower_segments[1].sequence == lower

        test_segments = benchmark(compactor.find_new_segments, test)
        display_segment_list(test_segments)

        assert len(test_segments) == 4
        assert test_segments[1].is_decision_kmer
        assert test_segments[1].sequence == test[:ksize]
        assert len(test_segments[1].sequence) == ksize
        assert len(test_segments[2].sequence) == len(test) - 1
        assert test_segments[2].left_anchor == graph.hash(core[pivot+1:pivot+ksize+1])
        assert test_segments[2].right_anchor == graph.hash(core[-ksize:])

    @using_ksize(15)
    @using_length(20)
    @pytest.mark.benchmark(group='cdbg-segments')
    def test_contained_loop(self, ksize, length, graph, compactor,
                                  circular, check_fp, benchmark):
        sequence = circular()
        check_fp()

        segments = benchmark(compactor.find_new_segments, sequence)
        assert len(segments) == 3
        assert segments[1].left_anchor == graph.hash(sequence[:ksize])
        assert segments[1].right_anchor == graph.hash(sequence[length-1:length+ksize-1])

    @using_ksize(15)
    @using_length(50)
    @pytest.mark.benchmark(group='cdbg-segments')
    def test_zero_new_segments(self, ksize, length, graph, compactor,
                                     linear_path, check_fp, benchmark):
        sequence = linear_path()
        check_fp()

        compactor.update_sequence(sequence)
        segments = benchmark(compactor.find_new_segments, sequence)
        assert len(segments) == 0


    @using_ksize(15)
    @using_length(100)
    @pytest.mark.benchmark(group='cdbg-segments')
    def test_two_or_more_decisions(self, ksize, length, graph, compactor,
                                           snp_bubble, check_fp, benchmark):
        (wild, mut), pivotL, pivotR = snp_bubble()
        print(wild)
        print(' ' * (pivotL - 1), wild[pivotL:pivotL+ksize])
        print(' ' * (pivotR - 1), wild[pivotR:pivotR+ksize])
        print(mut)

        wilds = [wild[:ksize], wild[pivotL+1:pivotR+ksize-1]]
        for wild in wilds:
            compactor.update_sequence(wild)
        assert compactor.cdbg.n_unodes == 2
        assert compactor.cdbg.n_dnodes == 0
        assert compactor.cdbg.n_unitig_ends == 3
        
        segments = benchmark(compactor.find_new_segments, mut)
        assert len(segments) == 7

        assert segments[1].start == 1
        assert segments[1].length == pivotL + ksize - 2
        assert segments[2].start == pivotL
        assert segments[2].length == ksize
        assert segments[3].start == pivotL + 1
        assert segments[3].length == (pivotR + ksize - 1) - (pivotL + 1)
        assert segments[4].start == pivotR
        assert segments[4].length == ksize
        assert segments[5].start == pivotR + 1
        assert segments[5].length == len(mut) - (pivotR + 1)



class TestDecisionNodes(object):

    @using_ksize(15)
    @using_length(100)
    def test_new_decision_from_fork(self, ksize, length, graph, compactor,
                                          left_fork, check_fp):
        '''New decision node of form (begin)-[D]-[S]-(end)
        '''

        (core, branch), pivot = left_fork()
        check_fp()

        upper = branch
        lower = core[:pivot+ksize-1]
        test = core[pivot:]

        compactor.update_sequence(upper)
        compactor.update_sequence(lower)
        assert compactor.cdbg.n_dnodes == 0

        compactor.update_sequence(test)
        assert compactor.cdbg.n_dnodes == 1
        dnode_hash = graph.hash(test[:ksize])
        print('dnode hash:', dnode_hash)
        assert compactor.cdbg.has_dnode(dnode_hash)

    @using_ksize(15)
    @using_length(100)
    def test_left_end_induced_decision_from_fork(self, ksize, length, graph, compactor,
                                             left_fork, check_fp):
        '''Decision node induced by segment end which is also end of sequence
           of form (begin)-[x]-[D]-[S]-(end)
        '''

        (core, branch), pivot = left_fork()
        check_fp()

        compactor.update_sequence(core)
        assert compactor.cdbg.n_dnodes == 0

        compactor.update_sequence(branch)
        assert compactor.cdbg.n_dnodes == 1
        dnode_hash = graph.hash(core[pivot:pivot+ksize])
        print('dnode hash:', dnode_hash)
        assert compactor.cdbg.has_dnode(dnode_hash)

    @using_ksize(15)
    @using_length(100)
    def test_right_end_induced_decision_from_fork(self, ksize, length, graph, compactor,
                                                        right_fork, check_fp):
        '''Decision node induced by segment end which is also end of sequence
           of form (begin)-[S]-[D]-[x]-(end)
        '''

        (core, branch), pivot = right_fork()
        check_fp()

        print(core, core[:pivot+1] + branch, sep='\n')
        print(core[pivot:pivot+ksize])

        compactor.update_sequence(core)
        assert compactor.cdbg.n_dnodes == 0

        compactor.update_sequence(branch)
        assert compactor.cdbg.n_dnodes == 1
        dnode_hash = graph.hash(core[pivot:pivot+ksize])
        print('dnode hash:', dnode_hash)
        assert compactor.cdbg.has_dnode(dnode_hash)

    @using_ksize(15)
    @using_length(100)
    def test_left_mid_induced_decision_from_fork(self, ksize, length, graph, compactor,
                                                       left_fork, check_fp):
        ''' Decision node is induced by a non-decision segment end,
            with flanking known sequence to its  right

            (begin)-[S]-[D x]-(end)
        '''
        (core, branch), pivot = left_fork()
        check_fp()

        branch = branch + core[pivot+ksize-1:]
        print(core, branch, sep='\n')
        print(core[pivot:pivot+ksize])

        compactor.update_sequence(core)
        assert compactor.cdbg.n_dnodes == 0

        compactor.update_sequence(branch)
        assert compactor.cdbg.n_dnodes == 1
        dnode_hash = graph.hash(core[pivot:pivot+ksize])
        print('dnode hash:', dnode_hash)
        assert compactor.cdbg.has_dnode(dnode_hash)

    @using_ksize(15)
    @using_length(100)
    def test_right_mid_induced_decision_from_fork(self, ksize, length, graph, compactor,
                                                        right_fork, check_fp):
        ''' Decision node is induced by a non-decision segment end,
            with flanking known sequence to its  right

            (begin)-[x D]-[S]-(end)
        '''
        (core, branch), pivot = right_fork()
        check_fp()

        branch = core[:pivot+1] + branch
        print(core, branch, sep='\n')
        print(core[pivot:pivot+ksize])

        compactor.update_sequence(core)
        assert compactor.cdbg.n_dnodes == 0

        compactor.update_sequence(branch)
        assert compactor.cdbg.n_dnodes == 1
        dnode_hash = graph.hash(core[pivot:pivot+ksize])
        print('dnode hash:', dnode_hash)
        assert compactor.cdbg.has_dnode(dnode_hash)

    @using_ksize(15)
    @using_length(100)
    def test_new_decision_node_segment_flanked(self, ksize, length, graph, compactor,
                                                     left_hairpin, check_fp):
        ''' Test flanked new decision node using a hairpin fixture, of form
            (begin)-[S]-[D]-[S]-[x]-(end) where [D] is the same k-mer as [x]
        '''
        sequence, pivot = left_hairpin()
        check_fp()
        print('d-kmer: ', sequence[pivot:pivot+ksize], graph.hash(sequence[pivot:pivot+ksize]))
        compactor.update_sequence(sequence)

        assert compactor.cdbg.n_dnodes == 1
        assert compactor.cdbg.has_dnode(graph.hash(sequence[pivot:pivot+ksize]))

    @using_ksize(15)
    @using_length(100)
    def test_trivial_unode_induced(self, ksize, length, graph, compactor,
                                         left_hairpin, check_fp):
        sequence, pivot = left_hairpin()
        check_fp()

        decision = sequence[pivot:pivot+ksize]
        compactor.update_sequence(decision)
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 1

        compactor.update_sequence(sequence)
        assert compactor.cdbg.has_dnode(graph.hash(decision))
        assert compactor.cdbg.query_unode_end(graph.hash(decision)) is None
        assert compactor.cdbg.n_unitig_ends == 4
        assert compactor.cdbg.n_unodes == 2


class TestUnitigBuildExtend(object):

    @using_ksize(15)
    @using_length(100)
    def test_left_fork_unode_creation(self, ksize, length, graph, compactor,
                                            left_fork, check_fp):
        '''New decision node of form (begin)-[D]-[S]-(end)
        '''

        (core, branch), pivot = left_fork()
        upper = branch
        lower = core[:pivot+ksize-1]
        test = core[pivot:]

        compactor.update_sequence(upper)
        assert compactor.cdbg.n_unodes == 1
        upper_unode = compactor.cdbg.query_unode_end(graph.hash(upper[:ksize]))
        assert upper_unode.sequence == upper

        compactor.update_sequence(lower)
        assert compactor.cdbg.n_dnodes == 0
        assert compactor.cdbg.n_unodes == 2
        lower_unode = compactor.cdbg.query_unode_end(graph.hash(lower[:ksize]))
        assert lower_unode.sequence == lower

        compactor.update_sequence(test)
        assert compactor.cdbg.n_dnodes == 1
        dnode_hash = graph.hash(test[:ksize])
        print('dnode hash:', dnode_hash)
        assert compactor.cdbg.has_dnode(dnode_hash)

        test_unode = compactor.cdbg.query_unode_end(graph.hash(test[1:ksize+1]))
        assert test_unode.sequence == test[1:]

    @using_ksize(15)
    @using_length(100)
    def test_extend_right(self, ksize, length, graph, compactor, linear_path, check_fp):
        sequence = linear_path()
        left = sequence[:length//2]

        compactor.update_sequence(left);
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.query_unode_end(graph.hash(left[:ksize])).sequence == left
        assert compactor.cdbg.n_unitig_ends == 2

        compactor.update_sequence(sequence)
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.query_unode_end(graph.hash(left[:ksize])).sequence == sequence
        assert compactor.cdbg.n_unitig_ends == 2

    @using_ksize(15)
    @using_length(100)
    def test_extend_left(self, ksize, length, graph, compactor, linear_path, check_fp):
        sequence = linear_path()
        check_fp()

        right = sequence[length//2:]

        compactor.update_sequence(right);
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.query_unode_end(graph.hash(right[:ksize])).sequence == right
        assert compactor.cdbg.n_unitig_ends == 2

        compactor.update_sequence(sequence)
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.query_unode_end(graph.hash(sequence[:ksize])).sequence == sequence
        assert compactor.cdbg.n_unitig_ends == 2

    @using_ksize(15)
    @using_length(100)
    def test_merge(self, ksize, length, graph, compactor, linear_path, check_fp):
        sequence = linear_path()
        check_fp()

        left = sequence[:length//2]
        right = sequence[length//2:]
        print(left)
        print(right)

        compactor.update_sequence(left);
        assert compactor.cdbg.n_unodes == 1
        unode = compactor.cdbg.query_unode_end(graph.hash(left[:ksize]))
        assert unode.sequence == left
        assert unode.meta == 'ISLAND'

        compactor.update_sequence(right)
        assert compactor.cdbg.n_unodes == 2
        unode = compactor.cdbg.query_unode_end(graph.hash(right[:ksize]))
        assert unode.sequence == right
        assert unode.meta == 'ISLAND'

        compactor.update_sequence(left + right)
        assert compactor.cdbg.n_unodes == 1
        unode = compactor.cdbg.query_unode_end(graph.hash(left[:ksize]))
        assert unode.sequence == left + right
        assert unode.meta == 'ISLAND'
        assert compactor.cdbg.query_unode_end(graph.hash(left[-ksize:])) is None
        assert compactor.cdbg.query_unode_end(graph.hash(right[:ksize])) is None

    @using_ksize(15)
    @using_length(100)
    def test_suffix_merge(self, ksize, length, graph, compactor, linear_path, check_fp):
        sequence = linear_path()
        check_fp()

        left = sequence[:length//2]
        right = sequence[length//2:]
        merger = left[-(ksize-1):] + right[:ksize-1]
        print(left)
        print(right)

        compactor.update_sequence(left);
        assert compactor.cdbg.n_unodes == 1
        unode = compactor.cdbg.query_unode_end(graph.hash(left[:ksize]))
        assert unode.sequence == left
        assert unode.meta == 'ISLAND'

        compactor.update_sequence(right)
        assert compactor.cdbg.n_unodes == 2
        unode = compactor.cdbg.query_unode_end(graph.hash(right[:ksize]))
        assert unode.sequence == right
        assert unode.meta == 'ISLAND'

        compactor.update_sequence(merger)
        assert compactor.cdbg.n_unodes == 1
        unode = compactor.cdbg.query_unode_end(graph.hash(left[:ksize]))
        assert unode.sequence == sequence
        assert unode.meta == 'ISLAND'
        assert compactor.cdbg.query_unode_end(graph.hash(left[-ksize:])) is None
        assert compactor.cdbg.query_unode_end(graph.hash(right[:ksize])) is None

    @using_ksize(7)
    @using_length(20)
    @pytest.mark.parametrize("offset", range(1,5))
    def test_overlap_merge(self, ksize, offset, length, graph, compactor, linear_path, check_fp):
        sequence = linear_path()
        check_fp()

        pivot = length // 2
        left = sequence[:pivot + offset]
        right = sequence[pivot:]
        print(left)
        print(right)

        compactor.update_sequence(left);
        assert compactor.cdbg.n_unodes == 1
        unode = compactor.cdbg.query_unode_end(graph.hash(left[:ksize]))
        assert unode.sequence == left
        assert unode.meta == 'ISLAND'

        compactor.update_sequence(right)
        assert compactor.cdbg.n_unodes == 2
        unode = compactor.cdbg.query_unode_end(graph.hash(right[:ksize]))
        assert unode.sequence == right
        assert unode.meta == 'ISLAND'

        compactor.update_sequence(sequence)
        assert compactor.cdbg.n_unodes == 1
        unode = compactor.cdbg.query_unode_end(graph.hash(left[:ksize]))
        assert unode.sequence == sequence
        assert unode.meta == 'ISLAND'
        assert compactor.cdbg.query_unode_end(graph.hash(left[-ksize:])) is None
        assert compactor.cdbg.query_unode_end(graph.hash(right[:ksize])) is None

    @using_ksize(15)
    @using_length(100)
    def test_trivial_merge_left(self, ksize, length, graph, compactor, linear_path, check_fp):
        sequence = linear_path()
        check_fp()

        left = sequence[:ksize]
        right = sequence[ksize:]
        merger = left[-(ksize-1):] + right[:ksize-1]
        print(left)
        print(right)

        compactor.update_sequence(left);
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 1
        unode = compactor.cdbg.query_unode_end(graph.hash(left[:ksize]))
        assert unode.sequence == left
        assert unode.meta == 'TRIVIAL'

        compactor.update_sequence(right)
        assert compactor.cdbg.n_unodes == 2
        assert compactor.cdbg.n_unitig_ends == 3
        unode = compactor.cdbg.query_unode_end(graph.hash(right[:ksize]))
        assert unode.sequence == right
        assert unode.meta == 'ISLAND'

        compactor.update_sequence(merger)
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 2
        unode = compactor.cdbg.query_unode_end(graph.hash(left[:ksize]))
        assert unode.sequence == sequence
        assert unode.meta == 'ISLAND'

    @using_ksize(15)
    @using_length(100)
    def test_trivial_merge_right(self, ksize, length, graph, compactor, linear_path, check_fp):
        sequence = linear_path()
        check_fp()

        left = sequence[:-ksize]
        right = sequence[-ksize:]
        merger = left[-(ksize-1):] + right[:ksize-1]
        print(left)
        print(right)

        compactor.update_sequence(left);
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 2
        unode = compactor.cdbg.query_unode_end(graph.hash(left[:ksize]))
        assert unode.sequence == left
        assert unode.meta == 'ISLAND'

        compactor.update_sequence(right)
        assert compactor.cdbg.n_unodes == 2
        assert compactor.cdbg.n_unitig_ends == 3
        unode = compactor.cdbg.query_unode_end(graph.hash(right[:ksize]))
        assert unode.sequence == right
        assert unode.meta == 'TRIVIAL'

        compactor.update_sequence(merger)
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 2
        unode = compactor.cdbg.query_unode_end(graph.hash(left[:ksize]))
        assert unode.sequence == sequence
        assert unode.meta == 'ISLAND'

    @using_length(50)
    @using_ksize(9)
    def test_suffix_extend(self, ksize, length, internal_pivot,
                                graph, compactor, linear_path, check_fp):
        sequence = linear_path()
        check_fp()

        pivot = internal_pivot
        left = sequence[:pivot+ksize]
        right = sequence[pivot+1:]
        print('\n', left, (' ' * (pivot + 1)) + right, sep='\n')

        compactor.update_sequence(left);
        assert compactor.cdbg.n_unodes == 1
        unode = compactor.cdbg.query_unode_end(graph.hash(left[:ksize]))
        assert unode.sequence == left
        assert unode.meta == 'ISLAND'

        compactor.update_sequence(right)
        assert compactor.cdbg.n_unodes == 1
        unode = compactor.cdbg.query_unode_end(graph.hash(right[-ksize:]))
        assert unode.sequence == sequence
        assert unode.meta == 'ISLAND'
        assert compactor.cdbg.query_unode_end(graph.hash(left[-ksize:])) is None


class TestUnitigSplit(object):

    @using_ksize(15)
    @using_length(50)
    def test_clip_from_left(self, right_sea, ksize, length, graph, compactor,
                                  check_fp):
        top, bottom = right_sea()
        check_fp()

        compactor.update_sequence(top)
        compactor.update_sequence(bottom)

        assert compactor.cdbg.n_dnodes == 1
        assert compactor.cdbg.n_unodes == 2
        assert compactor.cdbg.n_unitig_ends == 4

        assert compactor.cdbg.query_dnode(graph.hash(top[:ksize])).sequence == top[:ksize]
        assert compactor.cdbg.query_unode_end(graph.hash(top[:ksize])) is None
        
        top_unode = compactor.cdbg.query_unode_end(graph.hash(top[1:ksize+1]))
        assert top_unode is not None
        assert top_unode.sequence == top[1:]
        assert top_unode.right_end == graph.hash(top[-ksize:])
        assert top_unode.meta == 'TIP'
        
        bottom_unode = compactor.cdbg.query_unode_end(graph.hash(bottom[1:ksize+1]))
        assert bottom_unode is not None
        assert bottom_unode.sequence == bottom[1:]
        assert bottom_unode.right_end == graph.hash(bottom[-ksize:])
        assert bottom_unode.meta == 'TIP'

    @using_ksize(15)
    @using_length(50)
    def test_clip_from_right(self, left_sea, ksize, length, graph, compactor,
                                   check_fp):
        top, bottom = left_sea()
        check_fp()

        compactor.update_sequence(top)
        compactor.update_sequence(bottom)

        assert compactor.cdbg.n_dnodes == 1
        assert compactor.cdbg.n_unodes == 2
        assert compactor.cdbg.n_unitig_ends == 4

        assert compactor.cdbg.query_dnode(graph.hash(top[-ksize:])).sequence == top[-ksize:]
        assert compactor.cdbg.query_unode_end(graph.hash(top[-ksize:])) is None
        
        top_unode = compactor.cdbg.query_unode_end(graph.hash(top[-(ksize+1):-1]))
        assert top_unode is not None
        assert top_unode.sequence == top[:-1]
        assert top_unode.left_end == graph.hash(top[:ksize])
        assert top_unode.meta == 'TIP'
        
        bottom_unode = compactor.cdbg.query_unode_end(graph.hash(bottom[-(ksize+1):-1]))
        assert bottom_unode is not None
        assert bottom_unode.sequence == bottom[:-1]
        assert bottom_unode.left_end == graph.hash(bottom[:ksize])
        assert bottom_unode.meta == 'TIP'

    @using_ksize(15)
    @using_length(150)
    def test_induced_decision_to_unitig_extend(self, ksize, length, graph, compactor,
                                                     right_fork, check_fp):
        (core, branch), pivot = right_fork()
        check_fp()

        to_be_end_induced = core[ksize+2:pivot+1] + branch
        waist_left = core[:ksize+2]
        waist_right = core[pivot+ksize:]

        compactor.update_sequence(to_be_end_induced)
        assert compactor.cdbg.n_dnodes == 0
        assert compactor.cdbg.n_unodes == 1

        compactor.update_sequence(core)
        assert compactor.cdbg.n_dnodes == 1
        assert compactor.cdbg.n_unodes == 3
        assert compactor.cdbg.n_unitig_ends in [5,6]

        assert compactor.cdbg.query_dnode(graph.hash(core[pivot:pivot+ksize])).sequence == \
               core[pivot:pivot+ksize]

        assert compactor.cdbg.query_unode_end(graph.hash(core[:ksize])).right_end == \
               graph.hash(core[pivot-1:pivot+ksize-1])

        assert compactor.cdbg.query_unode_end(graph.hash(branch[:ksize])).right_end == \
               graph.hash(branch[-ksize:])
        
        assert compactor.cdbg.query_unode_end(graph.hash(core[pivot+1:pivot+ksize+1])).right_end == \
               graph.hash(core[-ksize:])
    
    @using_ksize(15)
    @using_length(100)
    def test_left_induced_split(self, ksize, length, graph, compactor,
                                      left_fork, check_fp):
        ''' Decision node is induced by a non-decision segment end,
            with flanking known sequence to its  right

            (begin)-[S]-[D x]-(end)
        '''
        (core, branch), pivot = left_fork()
        print(core, branch, sep='\n')
        print(core[pivot:pivot+ksize])

        compactor.update_sequence(core)
        assert compactor.cdbg.n_dnodes == 0
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 2

        compactor.update_sequence(branch)
        assert compactor.cdbg.n_dnodes == 1
        assert compactor.cdbg.n_unodes == 3
        assert compactor.cdbg.n_unodes == 3
        assert compactor.cdbg.n_unitig_ends in [5,6] # 5 with flank_left decision pivotition

        branch_unode = compactor.cdbg.query_unode_end(graph.hash(branch[:ksize]))
        assert branch_unode is not None
        assert branch_unode.sequence == branch
        assert branch_unode.left_end == graph.hash(branch[:ksize])
        assert branch_unode.right_end == graph.hash(branch[-ksize:])
    
        assert core[pivot:pivot+ksize] not in branch_unode.sequence
        assert compactor.cdbg.query_dnode(graph.hash(core[pivot:pivot+ksize])) is not None

        core_left_unode = compactor.cdbg.query_unode_end(graph.hash(core[:ksize]))
        assert core_left_unode is not None
        assert core_left_unode.sequence == core[:pivot+ksize-1]
        assert core_left_unode.left_end == graph.hash(core[:ksize])
        assert core_left_unode.right_end == graph.hash(core[pivot-1:pivot+ksize-1])

        core_right_unode = compactor.cdbg.query_unode_end(graph.hash(core[-ksize:]))
        assert core_right_unode.sequence == core[pivot+1:]
        assert core_right_unode.left_end == graph.hash(core[pivot+1:pivot+ksize+1])
        assert core_right_unode.right_end == graph.hash(core[-ksize:])


    @using_ksize(15)
    @using_length(100)
    def test_right_induced_split(self, ksize, length, graph, compactor,
                                       right_fork, check_fp):
        ''' Decision node is induced by a non-decision segment end,
            with flanking known sequence to its  right

            (begin)-[x D]-[S]-(end)
        '''
        (core, branch), pivot = right_fork()
        print(core, branch, sep='\n')
        print(core[pivot:pivot+ksize])

        compactor.update_sequence(core)
        assert compactor.cdbg.n_dnodes == 0
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 2

        compactor.update_sequence(branch)
        assert compactor.cdbg.n_dnodes == 1
        assert compactor.cdbg.n_unodes == 3
        assert compactor.cdbg.n_unitig_ends in [5,6]

        branch_unode = compactor.cdbg.query_unode_end(graph.hash(branch[:ksize]))
        assert branch_unode is not None
        assert branch_unode.sequence == branch
        assert branch_unode.left_end == graph.hash(branch[:ksize])
        assert branch_unode.right_end == graph.hash(branch[-ksize:])
    
        assert core[pivot:pivot+ksize] not in branch_unode.sequence
        assert compactor.cdbg.query_dnode(graph.hash(core[pivot:pivot+ksize])) is not None

        core_left_unode = compactor.cdbg.query_unode_end(graph.hash(core[:ksize]))
        assert core_left_unode is not None
        assert core_left_unode.sequence == core[:pivot+ksize-1]
        assert core_left_unode.left_end == graph.hash(core[:ksize])
        assert core_left_unode.right_end == graph.hash(core[pivot-1:pivot+ksize-1])

        core_right_unode = compactor.cdbg.query_unode_end(graph.hash(core[-ksize:]))
        assert core_right_unode.sequence == core[pivot+1:]
        assert core_right_unode.left_end == graph.hash(core[pivot+1:pivot+ksize+1])
        assert core_right_unode.right_end == graph.hash(core[-ksize:])

    @using_ksize(15)
    @using_length(100)
    def test_tandem_decision_unitig_clipping(self, ksize, length, graph, compactor,
                                          tandem_quad_forks, check_fp):
        (core, left_branches, right_branches), left_pivot, right_pivot = tandem_quad_forks()
        left_dkmer = core[left_pivot:left_pivot+ksize]
        right_dkmer = core[right_pivot:right_pivot+ksize]
        print('left d-node:', left_dkmer, left_pivot, graph.hash(left_dkmer))
        print('right d-node:', right_dkmer, right_pivot, graph.hash(right_dkmer))
        compactor.update_sequence(core)
        
        n_ends = 2
        n_unodes = 1
        for branch_num, branch in enumerate(left_branches):
            print('*** INSERT left branch', branch_num, file=sys.stderr)
            compactor.update_sequence(branch)
            if branch_num == 0:
                n_ends += 4 # the first branch induces the dnode and splits core
                n_unodes += 2
            else:
                n_ends += 2
                n_unodes += 1
            assert compactor.cdbg.n_dnodes == 1
            assert compactor.cdbg.n_unitig_ends == n_ends
            assert compactor.cdbg.n_unodes == n_unodes
            unode = compactor.cdbg.query_unode_end(graph.hash(branch[:ksize]))
            assert unode is not None
            assert unode.sequence == branch
            assert unode.right_end == graph.hash(branch[-ksize:])

        left_unode = compactor.cdbg.query_unode_end(graph.hash(core[:ksize]))
        assert left_unode is not None
        assert left_unode.right_end == graph.hash(core[left_pivot-1:left_pivot-1+ksize])
        assert left_unode.sequence == core[:left_pivot-1+ksize]

        for branch_num, branch in enumerate(right_branches):
            print('*** INSERT right branch', branch_num, file=sys.stderr)
            compactor.update_sequence(branch)
            n_ends += 2 # first branch induces the second d-node but clips
                        # the unitig rather than splitting it
            n_unodes += 1

            assert compactor.cdbg.n_dnodes == 2
            assert compactor.cdbg.n_unitig_ends == n_ends
            assert compactor.cdbg.n_unodes == n_unodes

            left_end = graph.hash(branch[:ksize])
            unode = compactor.cdbg.query_unode_end(left_end)
            print('left_end for branch is', left_end, file=sys.stderr)

            assert branch[:ksize] not in core
            branch_hashes = list(graph.hashes(branch[:ksize+3]))
            print(branch_hashes, file=sys.stderr)
            print("branch is length ", len(branch))

            assert unode is not None
            assert unode.sequence == branch
            assert unode.right_end == graph.hash(branch[-ksize:])

    @using_ksize(15)
    @using_length(100)
    def test_induced_chain(self, ksize, length, graph, compactor,
                                      snp_bubble, check_fp):

        (wild, snp), L, R = snp_bubble()
        check_fp()

        compactor.update_sequence(wild)
        assert compactor.cdbg.n_dnodes == 0
        assert compactor.cdbg.n_unodes == 1

        compactor.update_sequence(snp)
        assert compactor.cdbg.n_dnodes == 2
        assert compactor.cdbg.n_unodes == 4
        assert compactor.cdbg.n_unitig_ends == 8

        assert compactor.cdbg.query_dnode(graph.hash(wild[L:L+ksize])).sequence == \
               wild[L:L+ksize]
        assert compactor.cdbg.query_dnode(graph.hash(wild[R:R+ksize])).sequence == \
               wild[R:R+ksize]       
 
        left = compactor.cdbg.query_unode_end(graph.hash(wild[:ksize]))
        assert left is not None
        assert left.sequence == wild[:L+ksize-1]

        right = compactor.cdbg.query_unode_end(graph.hash(wild[-ksize:]))
        assert right is not None
        assert right.sequence == wild[R+1:]

        top = compactor.cdbg.query_unode_end(graph.hash(wild[L+1:L+1+ksize]))
        assert top is not None
        assert top.sequence == wild[L+1:R+ksize-1]

        bottom = compactor.cdbg.query_unode_end(graph.hash(snp[L+1:L+1+ksize]))
        assert bottom is not None
        assert bottom.sequence == snp[L+1:R+ksize-1]


    @using_ksize(15)
    @using_length(100)
    def test_induced_chain_hourglass(self, ksize, length, graph, compactor,
                                           hourglass_tangle, check_fp):

        (top, bottom), L = hourglass_tangle()
        check_fp()

        compactor.update_sequence(top)
        assert compactor.cdbg.n_dnodes == 0
        assert compactor.cdbg.n_unodes == 1

        compactor.update_sequence(bottom)
        assert compactor.cdbg.n_dnodes == 4
        assert compactor.cdbg.n_unodes == 4
        assert compactor.cdbg.n_unitig_ends == 8

        assert compactor.cdbg.query_dnode(graph.hash(top[L:L+ksize])) is not None
        assert compactor.cdbg.query_dnode(graph.hash(top[L+1:L+1+ksize])) is not None
        assert compactor.cdbg.query_dnode(graph.hash(bottom[L:L+ksize])) is not None
        assert compactor.cdbg.query_dnode(graph.hash(bottom[L+1:L+1+ksize])) is not None

        ltop = compactor.cdbg.query_unode_end(graph.hash(top[:ksize]))
        assert len(ltop) == L + ksize -1
        assert ltop.right_end == graph.hash(top[L-1:L-1+ksize])
        assert ltop.sequence == top[:L-1+ksize]

        rtop = compactor.cdbg.query_unode_end(graph.hash(top[-ksize:]))
        assert len(rtop) == length - L - 2
        assert rtop.left_end == graph.hash(top[L+2:L+2+ksize])
        assert rtop.sequence == top[L+2:]

        lbottom = compactor.cdbg.query_unode_end(graph.hash(bottom[:ksize]))
        assert len(lbottom) == L + ksize - 1
        assert lbottom.right_end == graph.hash(bottom[L-1:L-1+ksize])
        assert lbottom.sequence == bottom[:L-1+ksize]

        rbottom = compactor.cdbg.query_unode_end(graph.hash(bottom[-ksize:]))
        assert len(rbottom) == length - L - 2
        assert rbottom.left_end == graph.hash(bottom[L+2:L+2+ksize])
        assert rbottom.sequence == bottom[L+2:]

    @using_ksize(15)
    @using_length(100)
    def test_induced_bowtie_split(self, ksize, length, graph, compactor,
                                        bowtie_tangle, check_fp):

        (top, bottom), L = bowtie_tangle()
        check_fp()

        compactor.update_sequence(top)
        assert compactor.cdbg.n_dnodes == 0
        assert compactor.cdbg.n_unodes == 1

        compactor.update_sequence(bottom)
        assert compactor.cdbg.n_dnodes == 1
        assert compactor.cdbg.n_unodes == 4
        assert compactor.cdbg.n_unitig_ends == 8

        ltop = compactor.cdbg.query_unode_end(graph.hash(top[:ksize]))
        assert len(ltop) == L + ksize
        assert ltop.right_end == graph.hash(top[L:L+ksize])

        rtop = compactor.cdbg.query_unode_end(graph.hash(top[-ksize:]))
        assert len(rtop) == length - L - 2
        assert rtop.left_end == graph.hash(top[L+2:L+2+ksize])

        lbottom = compactor.cdbg.query_unode_end(graph.hash(bottom[:ksize]))
        assert len(lbottom) == L + ksize
        assert lbottom.right_end == graph.hash(bottom[L:L+ksize])

        rbottom = compactor.cdbg.query_unode_end(graph.hash(bottom[-ksize:]))
        assert len(rbottom) == length - L - 2
        assert rbottom.left_end == graph.hash(bottom[L+2:L+2+ksize])


class TestCircularUnitigs:

    @using_ksize(15)
    @using_length(20)
    def test_suffix_loop_at_end(self, ksize, length, graph, compactor,
                                      suffix_circular, check_fp):
        sequence = suffix_circular()
        check_fp()

        compactor.update_sequence(sequence)
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 1
        assert compactor.cdbg.n_dnodes == 0

        unode = compactor.cdbg.query_unode_end(graph.hash(sequence[:ksize]))
        assert unode is not None
        assert unode.right_end == graph.hash(sequence[:ksize])
        assert unode.meta == 'CIRCULAR'

    @using_ksize(15)
    @using_length(20)
    def test_contained_in_sequence(self, ksize, length, graph, compactor,
                                   circular, check_fp):
        sequence = circular()
        check_fp()

        compactor.update_sequence(sequence)
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 1
        unode = compactor.cdbg.query_unode_end(graph.hash(sequence[:ksize]))
        assert unode.right_end == graph.hash(sequence[:ksize])
        assert unode.sequence == sequence[:length+ksize-1]
        assert unode.meta == 'CIRCULAR'

    @using_ksize(15)
    @using_length(40)
    def test_circular_merge(self, ksize, length, graph, compactor,
                                                circular, check_fp):
        sequence = circular()
        check_fp()

        start = sequence[:length - 2]
        
        compactor.update_sequence(start)
        compactor.update_sequence(sequence)

        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 1

        unode = compactor.cdbg.query_unode_end(graph.hash(start[:ksize]))
        assert unode.right_end == graph.hash(start[:ksize])
        assert unode.sequence == sequence[:length+ksize-1]
        assert unode.meta == 'CIRCULAR'

    @using_ksize(7)
    @using_length(20)
    @pytest.mark.parametrize("offset", range(1,7), ids=lambda offset: 'offset={0}'.format(offset))
    def test_circular_cycling_merge(self, ksize, offset, length, graph, compactor,
                                          suffix_circular, check_fp):
        sequence = suffix_circular()
        check_fp()

        pivot = (length - ksize + 1) // 2
        start = sequence[:pivot]
        print(sequence, pivot)
        expected = sequence[pivot+offset:] + sequence[ksize-1:pivot+offset+ksize-1]
        
        compactor.update_sequence(start)
        compactor.update_sequence(sequence[pivot+offset:])

        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 2
        
        compactor.update_sequence(sequence)

        print('unitig should be:', expected)

        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 1

        unode = compactor.cdbg.query_unode_end(graph.hash(expected[:ksize]))
        assert unode.right_end == graph.hash(expected[:ksize])
        assert unode.sequence == expected
        assert unode.meta == 'CIRCULAR'

    @using_ksize(7)
    @using_length(20)
    @pytest.mark.parametrize("offset", range(1,6), ids=lambda offset: 'offset={0}'.format(offset))
    def test_circular_overlap_merge(self, ksize, offset, length, graph, compactor,
                                          suffix_circular, check_fp):
        sequence = suffix_circular()
        check_fp()

        start = sequence[:-ksize+1+offset]
        print(sequence, start)
        
        compactor.update_sequence(start)
        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 2

        compactor.update_sequence(sequence)

        assert compactor.cdbg.n_unodes == 1
        assert compactor.cdbg.n_unitig_ends == 1

        unode = compactor.cdbg.query_unode_end(graph.hash(sequence[:ksize]))
        assert unode.right_end == graph.hash(sequence[:ksize])
        assert unode.sequence == sequence
        assert unode.meta == 'CIRCULAR'

    @using_ksize(7)
    @using_length(21)
    def test_split_circular_tangle_chain(self, ksize, length, graph, compactor,
                                               suffix_circular_tangle, check_fp):
        (loop, inducer), lpivot = suffix_circular_tangle()
        check_fp()

        # create circular unitig
        # already tested
        compactor.update_sequence(loop)
        compactor.update_sequence(inducer)

        assert compactor.cdbg.n_unodes == 3
        assert compactor.cdbg.n_unitig_ends == 6
        assert compactor.cdbg.n_dnodes == 4

        assert compactor.cdbg.query_dnode(graph.hash(loop[lpivot:lpivot+ksize])).sequence == loop[lpivot:lpivot+ksize]
        assert compactor.cdbg.query_dnode(graph.hash(loop[lpivot+1:lpivot+1+ksize])).sequence == loop[lpivot+1:lpivot+1+ksize]
        
        cycled_unode = compactor.cdbg.query_unode_end(graph.hash(loop[lpivot-1:lpivot-1+ksize]))
        expected = loop[lpivot+2:] + loop[ksize-1:lpivot+ksize-1]
        assert cycled_unode is not None
        assert cycled_unode.sequence == expected
        assert cycled_unode.right_end == graph.hash(expected[-ksize:])

    @using_ksize(7)
    @using_length(20)
    @using_pivot(['flank_left', 'middle', 'flank_right'])
    def test_split_circular(self, ksize, length, graph, compactor,
                                  circular_key, check_fp):
        (loop, tail), pivot = circular_key()
        check_fp()

        compactor.update_sequence(loop)
        loop_unode = compactor.cdbg.query_unode_end(graph.hash(loop[:ksize]))
        assert loop_unode is not None
        assert loop_unode.right_end == graph.hash(loop[:ksize])
        assert loop_unode.sequence == loop
        assert loop_unode.meta == 'CIRCULAR'
        loop_unode = loop_unode.clone()

        compactor.update_sequence(tail)
        assert compactor.cdbg.n_dnodes == 1
        assert compactor.cdbg.n_unodes == 2
        assert compactor.cdbg.n_unitig_ends == 3

        cycled_loop_unode = compactor.cdbg.query_unode_end(graph.hash(loop[pivot-1:pivot-1+ksize]))
        assert cycled_loop_unode is not None
        assert cycled_loop_unode.left_end == graph.hash(loop[pivot+1:pivot+1+ksize])

        print('\n', loop_unode, sep='')
        print(cycled_loop_unode)

        print('\n', loop, sep='')
        print((' ' * pivot) + loop[pivot:pivot+ksize])
        print(loop_unode.sequence)
        print((' ' * (pivot+1)) + cycled_loop_unode.sequence)

    @using_ksize(7)
    @using_length(20)
    @using_pivot(['left', 'right'])
    def test_split_circular_on_end(self, ksize, length, graph, compactor,
                                         circular_key, check_fp):
        (loop, tail), pivot = circular_key()
        check_fp()
        print('length={0} pivot={1}'.format(length, pivot))
        print(loop, ' ' * pivot + tail)
        print(list(graph.hash(kmer) for kmer in kmers(loop, ksize)))

        compactor.update_sequence(loop)
        loop_unode = compactor.cdbg.query_unode_end(graph.hash(loop[:ksize]))
        assert loop_unode is not None
        assert loop_unode.right_end == graph.hash(loop[:ksize])
        assert loop_unode.sequence == loop
        assert loop_unode.meta == 'CIRCULAR'
        loop_unode = loop_unode.clone()

        compactor.update_sequence(tail)
        assert compactor.cdbg.n_dnodes == 1
        assert compactor.cdbg.n_unodes == 2
        assert compactor.cdbg.n_unitig_ends == 3

        dnode = compactor.cdbg.query_dnode(graph.hash(loop[pivot:pivot+ksize]))
        assert dnode.sequence == loop[pivot:pivot+ksize]

        if pivot == 0:
            # dnode is first k-mer in loop
            cycled_loop_unode = compactor.cdbg.query_unode_end(graph.hash(loop[pivot+1:pivot+1+ksize]))
            print(cycled_loop_unode)
            assert cycled_loop_unode is not None
            assert cycled_loop_unode.right_end == graph.hash(loop[-ksize:])
            assert cycled_loop_unode.meta == 'FULL'
        else:
            # dnode is last k-mer in loop
            cycled_loop_unode = compactor.cdbg.query_unode_end(graph.hash(loop[pivot-1:pivot-1+ksize]))   
            print(cycled_loop_unode)
            print('pivot:', pivot)
            assert cycled_loop_unode is not None
            assert cycled_loop_unode.left_end == graph.hash(loop[:ksize])
            assert cycled_loop_unode.meta == 'FULL'

        print('\n', loop_unode, sep='')
        print(cycled_loop_unode)

        print('\n', loop, sep='')
        print((' ' * pivot) + loop[pivot:pivot+ksize])
        print(loop_unode.sequence)
        print((' ' * (pivot+1)) + cycled_loop_unode.sequence)


class TestBreadthFirstTraversal:

    @using_ksize(21)
    @using_length(100)
    @pytest.mark.benchmark(group='cdbg-traversal')
    def test_single_component_from_left_unode(self, ksize, length, graph, compactor,
                                                    snp_bubble, check_fp, benchmark):

        (wild, snp), L, R = snp_bubble()
        check_fp()

        compactor.update_sequence(wild)
        compactor.update_sequence(snp)

        left_unode_root = compactor.cdbg.query_unode_end(graph.hash(wild[:ksize]))
        nodes = benchmark(compactor.cdbg.traverse_breadth_first, left_unode_root)
        assert len(nodes) == 6


    @using_ksize(21)
    @using_length(100)
    @pytest.mark.benchmark(group='cdbg-traversal')
    def test_single_component_from_right_unode(self, ksize, length, graph, compactor,
                                                     snp_bubble, check_fp, benchmark):

        (wild, snp), L, R = snp_bubble()
        check_fp()

        compactor.update_sequence(wild)
        compactor.update_sequence(snp)

        root = compactor.cdbg.query_unode_end(graph.hash(wild[-ksize:]))
        nodes = benchmark(compactor.cdbg.traverse_breadth_first, root)
        assert len(nodes) == 6

    @using_ksize(21)
    @using_length(100)
    @pytest.mark.benchmark(group='cdbg-traversal')
    def test_single_component_from_dnode(self, ksize, length, graph, compactor,
                                               snp_bubble, check_fp, benchmark):

        (wild, snp), L, R = snp_bubble()
        check_fp()

        compactor.update_sequence(wild)
        compactor.update_sequence(snp)

        root = compactor.cdbg.query_dnode(graph.hash(wild[L:L+ksize]))
        nodes = benchmark(compactor.cdbg.traverse_breadth_first, root)
        assert len(nodes) == 6


class TestFindConnectedComponents:

    @using_ksize(21)
    @using_length(100)
    @pytest.mark.benchmark(group='cdbg-find-components')
    def test_single_component(self, ksize, length, graph, compactor,
                                    snp_bubble, check_fp, benchmark):

        (wild, snp), L, R = snp_bubble()
        check_fp()

        compactor.update_sequence(wild)
        compactor.update_sequence(snp)

        components = benchmark(compactor.cdbg.find_connected_components)
        assert len(components) == 1
        print(components)
        assert len(components[0]) == 6
        dl, dr = graph.hash(wild[L:L+ksize]), graph.hash(wild[R:R+ksize])
        assert all((ID in components[0] for ID in (0,1,2,3,dl,dr)))

    @using_ksize(21)
    @using_length(100)
    @pytest.mark.parametrize('n_components', [10, 50, 100])
    @pytest.mark.benchmark(group='cdbg-find-components')
    def test_many_components(self, ksize, length, n_components, graph, compactor,
                                   snp_bubble, check_fp, benchmark):

        for _ in range(n_components):
            (wild, snp), L, R = snp_bubble()
            check_fp()

            compactor.update_sequence(wild)
            compactor.update_sequence(snp)

        components = benchmark(compactor.cdbg.find_connected_components)
        assert len(components) == n_components
