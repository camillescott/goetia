# boink/tests/test_usparsetagger.py
# Copyright (C) 2020 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest

from tests.utils import *

from boink import libboink
from boink.dbg import dBG
from boink.hashing import CanUnikmerShifter
from boink.storage import SparseppSetStorage


@pytest.fixture
def tagger(ksize):
    tagger_type = libboink.cdbg.USparseGraph[SparseppSetStorage].Graph
    storage = SparseppSetStorage.build()
    hasher = CanUnikmerShifter.build(ksize, 7)
    graph = dBG[type(storage), CanUnikmerShifter].build(storage, hasher)
    tagger = tagger_type.build(graph, hasher.ukhs_map)
    return tagger


@using(ksize=31, length=100)
class TestFindNewExtensions:

    def test_new_whole_sequence(self, ksize, length, tagger, linear_path):
        sequence = linear_path()
        extensions = tagger.find_new_extensions(sequence)

        assert len(extensions) == len(sequence) - ksize + 1
        for h, pos, (left, right) in extensions:
            assert len(left) == 4
            assert len(right) == 4