# boink/tests/test_cdbg.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
from boink.tests.utils import *


@pytest.fixture
def compactor(ksize, graph):
    from boink.cdbg import StreamingCompactor
    compactor = StreamingCompactor(graph)
    return compactor


@using_ksize(21)
@pytest.mark.parametrize('graph_type', ['BitStorage'], indirect=['graph_type'])
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
@pytest.mark.parametrize('graph_type', ['BitStorage'], indirect=['graph_type'])
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
@pytest.mark.parametrize('graph_type', ['BitStorage'], indirect=['graph_type'])
def test_find_decision_nodes_fork(ksize, graph, compactor, consumer, right_fork):
    (core, branch), pos = right_fork()
    print(core, ' ' * (pos + 1) + branch, sep='\n')

    positions, hashes = compactor.find_decision_nodes(core)
    assert positions == [pos]
    assert hashes == [graph.hash(core[pos:pos+ksize])]


@using_ksize(21)
@pytest.mark.parametrize('graph_type', ['BitStorage'], indirect=['graph_type'])
def test_find_decision_nodes_objects(ksize, graph, compactor, right_fork):
    (core, branch), pos = right_fork()
    print('\n', core, ' ' * (pos + 1) + branch, sep='\n')
    compactor.update(branch)
    compactor.update(core)

    dnode = list(compactor.get_cdbg_dnodes(core)).pop()
    assert dnode.in_degree == 1
    left_unode = list(compactor.cdbg.left_neighbors(dnode)).pop()
    print(left_unode.sequence)
    assert dnode.sequence == left_unode.sequence[-ksize:]

    assert dnode.out_degree == 2
    for unode in compactor.cdbg.right_neighbors(dnode):
        assert dnode.sequence == unode.sequence[:ksize]

    #assert node.out_degree == 2
    #assert hashes == [graph.hash(core[pos:pos+ksize])]
