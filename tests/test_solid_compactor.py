# boink/tests/test_cdbg.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import itertools
import sys

import pytest

from boink import libboink

from .utils import *
from .test_cdbg import compactor, compactor_type


@pytest.fixture(params=[1, 3])
def min_abund(request):
    return request.param


@pytest.fixture
def solid_compactor(graph, compactor, compactor_type, min_abund):
    
    _solid_compactor = compactor_type.SolidCompactor.build(compactor,
                                                           min_abund,
                                                           100000,
                                                           4)
    return _solid_compactor


@using(hasher_type=libboink.hashing.RollingHashShifter)
class TestFindSolidSegments:

    @using(ksize=21, length=100, min_abund=2)
    def test_no_segments(self, ksize, length, graph, compactor, solid_compactor,
                               min_abund, linear_path, check_fp):
        sequence = linear_path()

        segments = solid_compactor.find_solid_segments(sequence)
        assert len(segments) == 0

    @using(ksize=21, length=100)
    def test_one_segments(self, ksize, length, graph, compactor, solid_compactor,
                                min_abund, linear_path, check_fp):
        sequence = linear_path()

        for _ in range(min_abund):
            solid_compactor.insert_sequence(sequence)

        segments = solid_compactor.find_solid_segments(sequence)
        assert [tuple(segments[0])] == [(0, len(sequence))]
