# boink/tests/test_cdbg.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import itertools
import sys

import pytest

from boink.compactor import SolidStreamingCompactor
from boink.tests.utils import *
from boink.tests.test_cdbg import compactor


@pytest.fixture(params=[1, 3])
def min_abund(request):
    return request.param


@pytest.fixture
def solid_compactor(graph, compactor, min_abund):
    _solid_compactor = SolidStreamingCompactor.build(compactor,
                                                     min_abund,
                                                     graph.size,
                                                     graph.n_tables)
    return _solid_compactor


class TestFindSolidSegments:

    @using_ksize(7)
    @using_length(20)
    @pytest.mark.parametrize('min_abund', [2], indirect=True)
    def test_no_segments(self, ksize, length, graph, compactor, solid_compactor,
                               min_abund, linear_path, check_fp):
        sequence = linear_path()

        segments = solid_compactor.find_solid_segments(sequence)
        assert segments == []

    @using_ksize(7)
    @using_length(20)
    def test_one_segments(self, ksize, length, graph, compactor, solid_compactor,
                                min_abund, linear_path, check_fp):
        sequence = linear_path()

        for _ in range(min_abund):
            solid_compactor.update_sequence(sequence)

        segments = solid_compactor.find_solid_segments(sequence)
        assert segments == [(0, len(sequence))]
