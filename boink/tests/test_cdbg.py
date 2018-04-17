# boink/tests/test_cdbg.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import sys

import pytest
from boink.tests.utils import *


@pytest.fixture
def compactor(ksize, graph):
    from boink.cdbg import StreamingCompactor
    compactor = StreamingCompactor(graph)
    return compactor


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
    compactor.update(branch)
    compactor.update(core)

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
    compactor.update(top)
    compactor.update(bottom)
    compactor.update(core)

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
        compactor.update(top)
        compactor.update(bottom)
        if loop == 1:
            print('CORE update', file=sys.stderr)
            compactor.update(core)
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

    compactor.update(top1)
    compactor.update(bottom1)
    compactor.update(core1)

    dnode = list(compactor.get_cdbg_dnodes(core1)).pop()
    assert dnode.left_degree == 1
    left_unode = list(compactor.cdbg.left_neighbors(dnode)).pop()
    assert dnode.sequence == left_unode.sequence[-ksize:]

    assert dnode.right_degree == 3
    for unode in compactor.cdbg.right_neighbors(dnode):
        assert dnode.sequence == unode.sequence[:ksize]

    assert compactor.cdbg.n_dnodes == 1
    assert compactor.cdbg.n_unodes == 4

    compactor.update(top2)
    compactor.update(bottom2)
    print('CORE1 + CORE2 update', file=sys.stderr)
    compactor.update(core1 + core2)

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


@using_ksize(7)
@using_length(40)
@pytest.mark.parametrize('graph_type', ['_BitStorage'], indirect=['graph_type'])
def test_update_snp_bubble(ksize, graph, compactor, snp_bubble):
    (wildtype, snp), L, R = snp_bubble()

    compactor.update(wildtype)
    compactor.update(snp)
    
    assert compactor.cdbg.n_dnodes == 2
    assert compactor.cdbg.n_unodes == 4

    dnodes = list(compactor.get_cdbg_dnodes(wildtype))
    left, right = dnodes
    assert left.left_degree == 1
    assert left.right_degree == 2
    assert right.left_degree == 2
    assert right.right_degree == 1
