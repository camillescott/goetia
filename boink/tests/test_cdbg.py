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

from boink.cdbg import display_segment_list


@pytest.fixture
def compactor(ksize, graph):
    from boink.cdbg import StreamingCompactor
    compactor = StreamingCompactor(graph)
    return compactor


class TestFindNewSegments:

    @using_ksize(15)
    @using_length(100)
    @pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
    def test_fork_core_first(self, ksize, length, graph, compactor, right_fork,
                                   check_fp):
        (core, branch), pos = right_fork()
        print('INPUTS', core, core[:pos+1], ' ' * (pos + 1) + branch, sep='\n\n')

        core_segments = compactor.find_new_segments(core)
        branch_segments = compactor.find_new_segments(core[:pos+1] + branch)

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
        assert branch_segments[1].start == pos + 1
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
    @pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
    def test_right_decision_split(self, ksize, graph, compactor, right_fork,
                                  check_fp):
        (core, branch), pos = right_fork()
        print('INPUTS', core, core[:pos+1], ' ' * (pos + 1) + branch, sep='\n\n')

        branch_segments = compactor.find_new_segments(branch)
        core_segments = compactor.find_new_segments(core)

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
        assert graph.left_degree(core[pos:pos+ksize]) == 1
        assert graph.right_degree(core[pos:pos+ksize]) == 2

        assert core[:pos+ksize-1] == core_segments[1].sequence
        assert core[pos+1:] == core_segments[3].sequence


    @using_ksize(15)
    @using_length(100)
    @pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
    def test_merge_no_decisions(self, ksize, length, graph, compactor, linear_path, check_fp):
        sequence = linear_path()
        left = sequence[:length//2]
        right = sequence[length//2:]
        print(left)
        print(right)

        left_segments = compactor.find_new_segments(left)
        right_segments = compactor.find_new_segments(right)

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
        merged_segments = compactor.find_new_segments(merged)
        assert len(merged_segments) == 3
        assert merged_segments[0].is_null
        assert merged_segments[-1].is_null
        merged_segment = merged_segments[1]

        assert merged_segment.sequence == merged_new
        assert merged_segment.left_anchor == graph.hash(merged_new[:ksize])
        assert merged_segment.right_anchor == graph.hash(merged_new[-ksize:])


    @using_ksize(15)
    @using_length(100)
    @pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
    def test_right_decision_on_end(self, ksize, graph, compactor, right_fork,
                                         check_fp):
        (core, branch), pos = right_fork()
        # pos is start position of decision k-mer
        print('INPUTS', core, core[:pos+1], ' ' * (pos + 1) + branch, sep='\n\n')
        upper = branch
        lower = core[pos+1:]
        test = core[:pos+ksize]

        upper_segments = compactor.find_new_segments(upper)
        lower_segments = compactor.find_new_segments(lower)

        assert len(upper_segments) == 3
        assert upper_segments[1].sequence == upper
        assert len(lower_segments) == 3
        assert lower_segments[1].sequence == lower

        test_segments = compactor.find_new_segments(test)
        display_segment_list(test_segments)

        assert test_segments[0].is_null
        assert len(test_segments) == 4
        assert test_segments[2].is_decision_kmer
        assert len(test_segments[2].sequence) == ksize
        assert len(test_segments[1].sequence) == pos + ksize - 1
        assert test_segments[1].left_anchor == graph.hash(core[:ksize])
        assert test_segments[1].right_anchor == graph.hash(core[pos-1:pos+ksize-1])
        assert test_segments[-1].is_null


    @using_ksize(15)
    @using_length(100)
    @pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
    def test_left_decision_on_end(self, ksize, length, graph, compactor, left_fork,
                                        check_fp):
        (core, branch), pos = left_fork()
        # pos is start position of decision k-mer
        print('INPUTS', core, core[:pos+1], ' ' * (pos + 1) + branch, sep='\n\n')
        upper = branch
        lower = core[:pos+ksize-1]
        test = core[pos:]

        upper_segments = compactor.find_new_segments(upper)
        lower_segments = compactor.find_new_segments(lower)

        assert len(upper_segments) == 3
        assert upper_segments[1].sequence == upper
        assert len(lower_segments) == 3
        assert lower_segments[1].sequence == lower

        test_segments = compactor.find_new_segments(test)
        display_segment_list(test_segments)

        assert len(test_segments) == 4
        assert test_segments[1].is_decision_kmer
        assert test_segments[1].sequence == test[:ksize]
        assert len(test_segments[1].sequence) == ksize
        assert len(test_segments[2].sequence) == len(test) - 1
        assert test_segments[2].left_anchor == graph.hash(core[pos+1:pos+ksize+1])
        assert test_segments[2].right_anchor == graph.hash(core[-ksize:])


class TestDecisionNodes(object):

    @using_ksize(15)
    @using_length(100)
    @pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
    def test_new_decision_from_fork(self, ksize, length, graph, compactor,
                                          left_fork, check_fp):
        '''New decision node of form (begin)-[D]-[S]-(end)
        '''

        (core, branch), pos = left_fork()
        upper = branch
        lower = core[:pos+ksize-1]
        test = core[pos:]

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
    @pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
    def test_left_end_induced_decision_from_fork(self, ksize, length, graph, compactor,
                                             left_fork, check_fp):
        '''Decision node induced by segment end which is also end of sequence
           of form (begin)-[x]-[D]-[S]-(end)
        '''

        (core, branch), pos = left_fork()

        compactor.update_sequence(core)
        assert compactor.cdbg.n_dnodes == 0

        compactor.update_sequence(branch)
        assert compactor.cdbg.n_dnodes == 1
        dnode_hash = graph.hash(core[pos:pos+ksize])
        print('dnode hash:', dnode_hash)
        assert compactor.cdbg.has_dnode(dnode_hash)

    @using_ksize(15)
    @using_length(100)
    @pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
    def test_right_end_induced_decision_from_fork(self, ksize, length, graph, compactor,
                                                        right_fork, check_fp):
        '''Decision node induced by segment end which is also end of sequence
           of form (begin)-[S]-[D]-[x]-(end)
        '''

        (core, branch), pos = right_fork()
        print(core, core[:pos+1] + branch, sep='\n')
        print(core[pos:pos+ksize])

        compactor.update_sequence(core)
        assert compactor.cdbg.n_dnodes == 0

        compactor.update_sequence(branch)
        assert compactor.cdbg.n_dnodes == 1
        dnode_hash = graph.hash(core[pos:pos+ksize])
        print('dnode hash:', dnode_hash)
        assert compactor.cdbg.has_dnode(dnode_hash)

    @using_ksize(15)
    @using_length(100)
    @pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
    def test_left_mid_induced_decision_from_fork(self, ksize, length, graph, compactor,
                                                       left_fork, check_fp):
        ''' Decision node is induced by a non-decision segment end,
            with flanking known sequence to its  right

            (begin)-[S]-[D x]-(end)
        '''
        (core, branch), pos = left_fork()
        branch = branch + core[pos+ksize-1:]
        print(core, branch, sep='\n')
        print(core[pos:pos+ksize])

        compactor.update_sequence(core)
        assert compactor.cdbg.n_dnodes == 0

        compactor.update_sequence(branch)
        assert compactor.cdbg.n_dnodes == 1
        dnode_hash = graph.hash(core[pos:pos+ksize])
        print('dnode hash:', dnode_hash)
        assert compactor.cdbg.has_dnode(dnode_hash)

    @using_ksize(15)
    @using_length(100)
    @pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
    def test_right_mid_induced_decision_from_fork(self, ksize, length, graph, compactor,
                                                        right_fork, check_fp):
        ''' Decision node is induced by a non-decision segment end,
            with flanking known sequence to its  right

            (begin)-[x D]-[S]-(end)
        '''
        (core, branch), pos = right_fork()
        branch = core[:pos+1] + branch
        print(core, branch, sep='\n')
        print(core[pos:pos+ksize])

        compactor.update_sequence(core)
        assert compactor.cdbg.n_dnodes == 0

        compactor.update_sequence(branch)
        assert compactor.cdbg.n_dnodes == 1
        dnode_hash = graph.hash(core[pos:pos+ksize])
        print('dnode hash:', dnode_hash)
        assert compactor.cdbg.has_dnode(dnode_hash)



@using_ksize(21)
@pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
def test_insert_fork_noop(ksize, graph, compactor, right_fork):
    '''Nothing should happen with this ordering,
    because branch does not include the decision node and
    all the k-mers in core have already been seen when it's
    inserted a second time.'''

    (core, branch), pos = right_fork()
    print(core, ' ' * (pos + 1) + branch, sep='\n')

    positions, hashes = compactor.insert_sequence(core)
    assert positions == []
    assert hashes == []

    positions, hashes = compactor.insert_sequence(branch)
    assert positions == []
    assert hashes == []

    positions, hashes = compactor.insert_sequence(core)
    assert positions == []
    assert hashes == []


@using_ksize(21)
@pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
def test_insert_fork(ksize, graph, compactor, right_fork):
    (core, branch), pos = right_fork()
    print(core, ' ' * (pos + 1) + branch, sep='\n')

    positions, hashes = compactor.insert_sequence(branch)
    assert positions == []
    assert hashes == []

    positions, hashes = compactor.insert_sequence(core)
    assert positions == [pos]
    assert hashes == [graph.hash(core[pos:pos+ksize])]


@using_ksize(21)
@pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
def test_find_decision_kmers(ksize, graph, compactor, consumer, right_fork):
    (core, branch), pos = right_fork()
    print(core, ' ' * (pos + 1) + branch, sep='\n')

    positions, hashes = compactor.find_decision_kmers(core)
    assert positions == [pos]
    assert hashes == [graph.hash(core[pos:pos+ksize])]


@using_ksize(21)
@pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
def test_update_fork(ksize, graph, compactor, right_fork):
    (core, branch), pos = right_fork()
    print('\n', core, ' ' * (pos + 1) + branch, sep='\n')
    compactor.update_sequence(branch)
    compactor.wait_on_updates()
    assert compactor.cdbg.n_unodes == 1
    island = list(compactor.cdbg.unodes()).pop()
    assert island.meta == 'ISLAND'
    compactor.update_sequence(core)
    compactor.wait_on_updates()

    dnode = list(compactor.get_cdbg_dnodes(core)).pop()
    assert dnode.left_degree == 1
    left_unode = list(compactor.cdbg.left_neighbors(dnode)).pop()
    print(left_unode.sequence)
    assert dnode.sequence == left_unode.sequence[-ksize:]

    assert dnode.right_degree == 2
    for unode in compactor.cdbg.right_neighbors(dnode):
        assert dnode.sequence == unode.sequence[:ksize]

    assert compactor.cdbg.n_dnodes == 1
    assert compactor.cdbg.n_unodes == 3


@using_ksize(21)
@pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
def test_update_triple_fork(ksize, graph, compactor, right_triple_fork):
    (core, top, bottom), pos = right_triple_fork()
    print('\n', core, ' ' * (pos + 1) + top, sep='\n')
    compactor.update_sequence(top)
    compactor.update_sequence(bottom)
    compactor.update_sequence(core)
    compactor.wait_on_updates()

    dnode = list(compactor.get_cdbg_dnodes(core)).pop()
    assert dnode.left_degree == 1
    left_unode = list(compactor.cdbg.left_neighbors(dnode)).pop()
    print(left_unode.sequence)
    assert dnode.sequence == left_unode.sequence[-ksize:]

    assert dnode.right_degree == 3
    for unode in compactor.cdbg.right_neighbors(dnode):
        assert dnode.sequence == unode.sequence[:ksize]

    assert compactor.cdbg.n_dnodes == 1
    assert compactor.cdbg.n_unodes == 4


@using_ksize(21)
@pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
def test_double_update(ksize, graph, compactor, right_triple_fork,
                        random_sequence):
    (core, top, bottom), pos = right_triple_fork()
    print('\n', core, ' ' * (pos + 1) + top, sep='\n')


    for loop in (1,2):
        compactor.update_sequence(top)
        compactor.update_sequence(bottom)
        compactor.wait_on_updates()
        if loop == 1:
            print('CORE update', file=sys.stderr)
            compactor.update_sequence(core)
        if loop == 2:
            print('CORE2 update', file=sys.stderr)
            core = random_sequence() + core

        dnode = list(compactor.get_cdbg_dnodes(core)).pop()
        assert dnode.left_degree == 1
        left_unode = list(compactor.cdbg.left_neighbors(dnode)).pop()
        print(left_unode.sequence)
        assert dnode.sequence == left_unode.sequence[-ksize:]

        assert dnode.right_degree == 3
        for unode in compactor.cdbg.right_neighbors(dnode):
            assert dnode.sequence == unode.sequence[:ksize]

        assert compactor.cdbg.n_dnodes == 1
        assert compactor.cdbg.n_unodes == 4


@using_ksize(21)
@pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
def test_component_merge(ksize, graph, compactor, right_triple_fork,
                        random_sequence):
    (core1, top1, bottom1), pos = right_triple_fork()
    (core2, top2, bottom2), pos = right_triple_fork()

    compactor.update_sequence(top1)
    compactor.update_sequence(bottom1)
    compactor.update_sequence(core1)
    compactor.wait_on_updates()

    dnode = list(compactor.get_cdbg_dnodes(core1)).pop()
    assert dnode.left_degree == 1
    left_unode = list(compactor.cdbg.left_neighbors(dnode)).pop()
    assert dnode.sequence == left_unode.sequence[-ksize:]

    assert dnode.right_degree == 3
    for unode in compactor.cdbg.right_neighbors(dnode):
        assert dnode.sequence == unode.sequence[:ksize]

    assert compactor.cdbg.n_dnodes == 1
    assert compactor.cdbg.n_unodes == 4

    compactor.update_sequence(top2)
    compactor.update_sequence(bottom2)
    print('CORE1 + CORE2 update', file=sys.stderr)
    compactor.update_sequence(core1 + core2)
    compactor.wait_on_updates()

    assert compactor.cdbg.n_dnodes == 2
    assert compactor.cdbg.n_unodes == 7

    dnodes = list(compactor.get_cdbg_dnodes(core1 + core2))
    for dnode in dnodes:
        assert dnode.left_degree == 1
        assert dnode.right_degree == 3

    right  = dnodes[1]
    middle = list(compactor.cdbg.left_neighbors(right)).pop()
    left   = list(compactor.cdbg.left_neighbors(middle)).pop()

    assert left == dnodes[0]
    assert left.sequence == middle.sequence[:ksize]
    assert middle.sequence[-ksize:] == right.sequence


@using_ksize(21)
@pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
def test_update_snp_bubble(ksize, graph, compactor, snp_bubble):
    (wildtype, snp), L, R = snp_bubble()

    compactor.update_sequence(wildtype)
    compactor.update_sequence(snp)
    compactor.wait_on_updates()
    
    assert compactor.cdbg.n_dnodes == 2
    assert compactor.cdbg.n_unodes == 4

    dnodes = list(compactor.get_cdbg_dnodes(wildtype))
    left, right = dnodes
    assert left.left_degree == 1
    assert left.right_degree == 2
    assert right.left_degree == 2
    assert right.right_degree == 1


@using_ksize(21)
@pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
def test_update_trivial_dnodes(ksize, graph, compactor, tandem_quad_forks):
    (core, left_branches, right_branches), S_l, S_r = tandem_quad_forks()
    
    for branch in itertools.chain(left_branches, right_branches):
        compactor.update_sequence(branch)
    print('CORE update', file=sys.stderr)
    compactor.update_sequence(core)
    compactor.wait_on_updates()

    dnodes = list(compactor.get_cdbg_dnodes(core))
    left, right = dnodes
    assert left.left_degree == 1
    assert left.right_degree == 4
    assert right.left_degree == 1
    assert right.right_degree == 4
    assert left.sequence[1:] == right.sequence[:-1]

    trivial_unode = list(compactor.cdbg.left_neighbors(right)).pop()
    assert trivial_unode.sequence[:-1] == left.sequence
    assert trivial_unode.sequence[1:] == right.sequence
    assert len(trivial_unode) == ksize + 1


@using_ksize(21)
@pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
def test_linear_merge(ksize, graph, compactor, linear_path):
    seq1 = linear_path()
    seq2 = linear_path()

    compactor.update_sequence(seq1)
    compactor.wait_on_updates()
    unode = list(compactor.cdbg.unodes()).pop()
    assert len(unode) == len(seq1)
    assert unode.sequence == seq1
    assert compactor.cdbg.n_unodes == 1

    compactor.update_sequence(seq2)
    compactor.wait_on_updates()
    assert compactor.cdbg.n_unodes == 2

    unode1, unode2 = list(compactor.cdbg.unodes())
    if unode1.node_id > unode2.node_id:
        unode2, unode1 = unode1, unode2
    assert unode1.sequence == seq1
    assert unode2.sequence == seq2

    compactor.update_sequence(seq1 + seq2)
    compactor.wait_on_updates()
    assert compactor.cdbg.n_unodes == 1
    unode = list(compactor.cdbg.unodes()).pop()
    assert unode.sequence == seq1 + seq2


@using_ksize(9)
@using_length(50)
@pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
def test_tip_extend(ksize, graph, compactor, consumer, right_fork, linear_path):
    (core, branch), pos = right_fork()
    compactor.update_cdbg(branch)
    compactor.update_cdbg(core)
    compactor.wait_on_updates()
    assert compactor.cdbg.n_dnodes == 1
    assert compactor.cdbg.n_unodes == 3

    seq = linear_path()
    print('Add ISLAND', seq, file=sys.stderr)
    compactor.update_cdbg(seq)
    compactor.wait_on_updates()
    assert compactor.cdbg.n_dnodes == 1
    assert compactor.cdbg.n_unodes == 4

    print('MERGE island with', seq+core, file=sys.stderr)
    compactor.update_sequence(seq + core)
    compactor.wait_on_updates()
    assert compactor.cdbg.n_dnodes == 1
    assert compactor.cdbg.n_unodes == 3
    
    dnode = next(compactor.cdbg.dnodes())
    extended = next(compactor.cdbg.left_neighbors(dnode))
    assert extended.sequence == seq + core[:pos+ksize]


@using_ksize(21)
@using_length(70)
@pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
def test_multi_update(ksize, graph, compactor, consumer,
                      right_fork, linear_path, snp_bubble):
    (wildtype, snp), L, R = snp_bubble()

    compactor.update_cdbg(wildtype)
    compactor.update_cdbg(snp)
    compactor.wait_on_updates()
    
    assert compactor.cdbg.n_dnodes == 2
    assert compactor.cdbg.n_unodes == 4

    dnodes = list(compactor.get_cdbg_dnodes(wildtype))
    left, right = dnodes
    assert left.left_degree == 1
    assert left.right_degree == 2
    assert right.left_degree == 2
    assert right.right_degree == 1

    wildtype_extend = wildtype[-ksize:] + linear_path()

    compactor.update_cdbg(wildtype_extend)
    compactor.wait_on_updates()

    assert compactor.cdbg.n_dnodes == 2
    assert compactor.cdbg.n_unodes == 4

    dnodes = list(compactor.get_cdbg_dnodes(wildtype_extend))
    assert len(dnodes) == 0

    (core, branch), pos = right_fork()
    compactor.update_cdbg(branch)
    compactor.wait_on_updates()
    assert compactor.cdbg.n_unodes == 7

    compactor.update_cdbg(core)
    compactor.wait_on_updates()
    dnode = list(compactor.get_cdbg_dnodes(core)).pop()
    assert dnode.left_degree == 1
    assert dnode.right_degree == 2
    assert compactor.cdbg.n_dnodes == 3

    compactor.update_sequence(wildtype_extend + core)
    compactor.wait_on_updates()

    assert compactor.cdbg.n_dnodes == 3
    assert compactor.cdbg.n_unodes == 6
