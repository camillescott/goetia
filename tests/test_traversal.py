#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : test_traversal.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 08.10.2019


import pytest

from boink import libboink
from boink.traversal import Assembler, STATES

from .utils import *


@pytest.fixture
def asm(request, ksize, graph):

    return Assembler(graph)



def test_assembler_cursor(asm, ksize):
    seed = 'A' * ksize

    asm.cursor = seed
    assert asm.cursor == seed


def test_assembler_cursor_wrong_size(asm, ksize):
    seed = 'A' * (ksize - 1)

    with pytest.raises(TypeError):
        asm.cursor = seed


class TestLinear:

    def test_all_start_positions(self, ksize, linear_path, asm, consume, check_fp):
        # assemble entire contig, starting from wherever
        contig = linear_path()
        check_fp()
        consume()

        for start in range(0, len(contig), 150):
            if len(contig) - start < ksize:
                continue
            asm.clear_seen()
            path, (lstate, lend), (rstate, rend) = asm.assemble(contig[start:start + ksize])
            assert path == contig, (len(path), len(contig), start)
            assert lstate == rstate == STATES.STOP_FWD

    @using_ksize(9)
    @using_length(81)
    def test_all_left_to_beginning(self, ksize, length, linear_path, asm, consume, check_fp):
        # assemble directed left
        contig = linear_path()
        check_fp()
        consume()

        for start in range(0, len(contig), length // 5):
            if len(contig) - start < ksize:
                continue
            asm.clear_seen()
            path, (state, end) = asm.assemble_left(contig[start:start + ksize])
            assert path == contig[:start + ksize], start
            assert state == STATES.STOP_FWD

    @using_ksize(9)
    @using_length(81)
    def test_all_right_to_end(self, ksize, length, linear_path, asm, consume, check_fp):
        # assemble directed right
        contig = linear_path()
        check_fp()
        consume()

        for start in range(0, len(contig), length // 5):
            if len(contig) - start < ksize:
                continue
            asm.clear_seen()
            path, (state, end) = asm.assemble_right(contig[start:start + ksize])
            assert path == contig[start:], start
            assert state == STATES.STOP_FWD

    def test_circular(self, ksize, circular, asm, consume, check_fp):
        contig = circular()
        check_fp()
        consume()

        path, (state, end) = asm.assemble_right(contig[:ksize])
        print(path, ',', contig)
        assert path == contig[:len(path)]
        assert state == STATES.STOP_SEEN


class TestDecisions:

    def test_decision_fwd(self, ksize, right_fork, asm, graph, consume, check_fp):
        (sequence, branch), S = right_fork()
        check_fp()
        consume()

        path, (state, end) = asm.assemble_right(sequence[:ksize])
        assert state == STATES.DECISION_FWD
        assert path == sequence[:S+ksize]
        assert end == graph.hash(sequence[S:S+ksize])
        assert graph.right_degree(path[-ksize:]) == 2

        assembled_branch, (lstate, lend), (rstate, rend) = asm.assemble(branch[-ksize:])
        assert branch == assembled_branch

    
    def test_decision_rc(self, ksize, right_fork, asm, graph, consume, check_fp):
        '''Test that we assemble through fork from the right when we haven't
        assembled the core path already.'''
        (sequence, branch), S = right_fork()
        check_fp()
        consume()
        
        path, (lstate, lend), (rstate, rend) = asm.assemble(branch[-ksize:])

        assert branch == path
        assert graph.left_degree(path[:ksize]) == 1
        assert graph.right_degree(path[:ksize]) == 1
        assert lstate == STATES.DECISION_RC
        assert lend == graph.hash(branch[:ksize])
        assert rstate == STATES.STOP_FWD
        assert rend == graph.hash(branch[-ksize:])


    def test_triple_decision_fwd(self, ksize, right_triple_fork,
                                 asm, consume, check_fp):
        (core, top, bottom), S = right_triple_fork()
        check_fp()
        consume()
        
        path, (state, end) = asm.assemble_right(core[:ksize])
        assert path == core[:S+ksize]
        
        path, (lstate, lend), (rstate, rend) = asm.assemble(top[-ksize:])
        assert top == path

        path, (lstate, lend), (rstate, rend) = asm.assemble(bottom[-ksize:])
        assert bottom == path
