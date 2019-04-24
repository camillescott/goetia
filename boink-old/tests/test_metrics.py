# boink/tests/test_metrics.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest

from boink.metrics import ReservoirSample_Int as ReservoirSample

class TestReservoirSample:

    sample_sizes = [5, 10, 100, 1000]
    
    @pytest.mark.parametrize('sample_size', sample_sizes)
    def test_sample_size(self, sample_size):
        S = ReservoirSample(sample_size)
        assert S.sample_size == sample_size

    def test_n_sampled(self):
        S = ReservoirSample(10)
        for i in range(0, 100):
            assert S.n_sampled == i
            S.sample(7)

    @pytest.mark.parametrize('sample_size', sample_sizes)
    def test_result(self, sample_size):
        S = ReservoirSample(sample_size)
        for i in range(sample_size * 2):
            S.sample(7)
        assert S.n_sampled == sample_size * 2
        assert all((r == 7 for r in S.get_result()))
        assert len(S.get_result()) == sample_size
        

    @pytest.mark.parametrize('sample_size', sample_sizes)
    def test_clear(self, sample_size):
        S = ReservoirSample(sample_size)
        for i in range(sample_size * 2):
            S.sample(7)
        assert S.n_sampled == sample_size * 2
        assert all((r == 7 for r in S.get_result()))
        assert len(S.get_result()) == sample_size
        
        S.clear()
        assert S.n_sampled == 0
        assert all((r == 0 for r in S.get_result()))
        
    @pytest.mark.parametrize('sample_size', sample_sizes)
    def test_result_fewer_than_size(self, sample_size):
        S = ReservoirSample(sample_size)  
        for i in range(sample_size // 2):
            S.sample(7)
        counts = {}
        for r in S.get_result():
            counts[r] = counts.get(r,0) + 1
        assert counts[7] == sample_size // 2
        assert counts[0] == sample_size - (sample_size // 2)
