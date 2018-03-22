# boink/tests/test_assembly.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest

from boink.dbg import make_dBG
from boink.assembly import make_assembler
from boink.tests.utils import *
from boink.tests.test_dbg import dbg_type
from khmer.tests.graph_structure_fixtures import *


@pytest.fixture
def asm_type(request, ksize, dbg_type):

    def build():
        G = dbg_type()
        return make_assembler(G)

    return build


def test_assembler_type(asm_type):
    asm = asm_type()

    assert asm.storage_type == asm.Graph.storage_type
    assert asm.shifter_type == asm.Graph.shifter_type

    _, _, asm_suffix = type(asm).__name__.partition('_')
    _, _, dbg_suffix = type(asm.Graph).__name__.partition('_')

    assert asm_suffix == dbg_suffix, (asm_suffix, dbg_suffix)


def test_assembler_cursor(asm_type, ksize):
    seed = 'A' * ksize
    asm = asm_type()

    asm.cursor = seed
    assert asm.cursor == seed


def test_assembler_cursor_wrong_size(asm_type, ksize):
    seed = 'A' * (ksize - 1)
    asm = asm_type()

    with pytest.raises(ValueError):
        asm.cursor = seed


class TestNonBranching:

    def test_all_start_positions(self, ksize, linear_structure, asm_type):
        # assemble entire contig, starting from wherever
        _, contig = linear_structure()
        asm = asm_type()
        asm.Graph.add_sequence(contig)

        for start in range(0, len(contig), 150):
            asm.clear_seen()
            path = asm.assemble(contig[start:start + ksize])
            assert path == contig, (len(path), len(contig), start)

    def test_all_left_to_beginning(self, ksize, linear_structure, asm_type):
        # assemble directed left
        _, contig = linear_structure()
        asm = asm_type()
        asm.Graph.add_sequence(contig)


        for start in range(0, len(contig), 150):
            asm.clear_seen()
            path = asm.assemble_left(contig[start:start + ksize])
            print(path, ', ', contig[:start])
            assert path == contig[:start + ksize], start

    def test_all_right_to_end(self, ksize, linear_structure, asm_type):
        # assemble directed right
        _, contig = linear_structure()
        asm = asm_type()
        asm.Graph.add_sequence(contig)

        for start in range(0, len(contig), 150):
            asm.clear_seen()
            path = asm.assemble_right(contig[start:start + ksize])
            print(path, ', ', contig[:start])
            assert path == contig[start:], start

    def test_circular(self, ksize, circular_linear_structure, asm_type):
        _, contig = circular_linear_structure()
        asm = asm_type()
        asm.Graph.add_sequence(contig)

        path = asm.assemble_right(contig[:ksize])
        print(path, ',', contig)
        assert path == contig[:len(path)]



