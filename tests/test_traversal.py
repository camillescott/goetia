#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : test_traversal.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 08.10.2019


import pytest

from boink import libboink
from boink.traversal import STATES

from .utils import *


def test_cursor(graph, ksize):
    seed = 'A' * ksize

    graph.cursor = seed
    assert graph.cursor == seed


def test_cursor_wrong_size(graph, ksize):
    seed = 'A' * (ksize - 1)
    
    with pytest.raises(BaseException):
        graph.cursor = seed


class TestWalk:

    def test_to_string(self, ksize, linear_path, graph, consume, check_fp):
        contig = linear_path()
        check_fp()
        consume()

        walk = graph.walk_right(contig[0:ksize])
        assert walk.to_string() == contig
        assert str(walk) == contig

    def test_tail(self, ksize, linear_path, graph, consume, check_fp):
        contig = linear_path()
        check_fp()
        consume()

        walk = graph.walk_right(contig[0:ksize])
        assert walk.tail() == graph.get_hash()

    def test_glue(self, ksize, linear_path, graph, consume, check_fp):
        contig = linear_path()
        check_fp()
        consume()

        lwalk, rwalk = graph.walk(contig[1:ksize+1])
        assert lwalk.glue(rwalk) == contig
        assert rwalk.glue(lwalk) == contig


class TestLinear:

    def test_all_start_positions(self, ksize, linear_path, graph, consume, check_fp):
        # assemble entire contig, starting from wherever
        contig = linear_path()
        check_fp()
        consume()

        for start in range(0, len(contig), 150):
            if len(contig) - start < ksize:
                continue
            graph.clear_seen()
            lwalk, rwalk = graph.walk(contig[start:start + ksize])
            path = lwalk.glue(rwalk)

            assert path == contig, (len(path), len(contig), start)
            assert lwalk.end_state == rwalk.end_state == STATES.STOP_FWD

    @using(ksize=21, length=81)
    def test_all_left_to_beginning(self, ksize, length, linear_path, graph, consume, check_fp):
        # assemble directed left
        contig = linear_path()
        check_fp()
        consume()

        for start in range(0, len(contig), length // 5):
            if len(contig) - start < ksize:
                continue
            graph.clear_seen()
            walk = graph.walk_left(contig[start:start + ksize])
            assert walk.to_string() == contig[:start + ksize], start
            assert walk.end_state == STATES.STOP_FWD

    @using(ksize=21, length=81)
    def test_all_right_to_end(self, ksize, length, linear_path, graph, consume, check_fp):
        # assemble directed right
        contig = linear_path()
        check_fp()
        consume()

        for start in range(0, len(contig), length // 5):
            if len(contig) - start < ksize:
                continue
            graph.clear_seen()
            walk = graph.walk_right(contig[start:start + ksize])
            assert walk.to_string() == contig[start:], start
            assert walk.end_state == STATES.STOP_FWD

    def test_circular(self, ksize, circular, graph, consume, check_fp):
        contig = circular()
        check_fp()
        consume()

        walk = graph.walk_right(contig[:ksize])
        path = walk.to_string()
        print(path, ',', contig)
        assert path == contig[:len(path)]
        assert walk.end_state == STATES.STOP_SEEN


class TestDecisions:


    def test_decision_fwd(self, ksize, right_fork, graph, consume, check_fp):
        (sequence, branch), S = right_fork()
        check_fp()
        consume()

        walk = graph.walk_right(sequence[:ksize])
        path = walk.to_string()
        assert walk.end_state == STATES.DECISION_FWD
        assert path == sequence[:S+ksize]
        assert walk.tail() == graph.hash(sequence[S:S+ksize])
        assert graph.right_degree(path[-ksize:]) == 2

        lwalk, rwalk = graph.walk(branch[-ksize:])
        assert branch == lwalk.glue(rwalk)

    
    def test_decision_rev(self, ksize, right_fork, graph, consume, check_fp):
        '''Test that we assemble through fork from the right when we haven't
        assembled the core path already.'''
        (sequence, branch), S = right_fork()
        check_fp()
        consume()
        
        lwalk, rwalk = graph.walk(branch[-ksize:])
        path = lwalk.glue(rwalk)

        assert branch == path
        assert graph.left_degree(path[:ksize]) == 1
        assert graph.right_degree(path[:ksize]) == 1
        assert lwalk.end_state == STATES.DECISION_BKW
        print('Hash:', graph.hash(sequence[S:S+ksize]))
        assert lwalk.tail() == graph.hash(branch[:ksize])
        assert rwalk.end_state == STATES.STOP_FWD
        assert rwalk.tail() == graph.hash(branch[-ksize:])

    def test_start_from_reverse_decision_left(self, ksize, right_fork, graph, consume, check_fp):
        (sequence, branch), S = right_fork()
        check_fp()
        consume()

        walk = graph.walk_left(sequence[S:S+ksize])
        assert walk.to_string() == sequence[:S+ksize]
        assert walk.end_state == STATES.STOP_FWD

    def test_start_from_reverse_decision_right(self, ksize, left_fork, graph, consume, check_fp):
        (sequence, branch), S = left_fork()
        check_fp()
        consume()

        walk = graph.walk_right(sequence[S:S+ksize])
        assert walk.to_string() == sequence[S:]
        assert walk.end_state == STATES.STOP_FWD

    def test_triple_decision_fwd(self, ksize, right_triple_fork,
                                 graph, consume, check_fp):
        (core, top, bottom), S = right_triple_fork()
        check_fp()
        consume()
        
        rwalk = graph.walk_right(core[:ksize])
        assert str(rwalk) == core[:S+ksize]
        
        lwalk, rwalk = graph.walk(top[-ksize:])
        assert top == lwalk.glue(rwalk)

        lwalk, rwalk = graph.walk(bottom[-ksize:])
        assert bottom == lwalk.glue(rwalk)
