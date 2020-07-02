# goetia/tests/test_saturation.py
# Copyright (C) 2020 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from pprint import pprint
from statistics import mean

import numpy as np
import pytest

from goetia.saturation import (SlidingWindow, SlidingCutoff, normalized_mean,
                               all_cutoff, median_cutoff)


class TestSlidingWindow:

    def test_mean_without_time(self):
        vals = [1, 2, 3, 4, 5]
        window = SlidingWindow(3, mean)

        result, time = window.push(vals[0])
        assert np.isnan(result)
        assert time == 0

        result, time = window.push(vals[1])
        assert np.isnan(result)
        assert time == 1

        result, time = window.push(vals[2])
        assert result == 2
        assert time == 2

        result, time = window.push(vals[3])
        assert result == 3
        assert time == 3

        result, time = window.push(vals[4])
        assert result == 4
        assert time == 4

    def test_mean_with_time(self):
        vals = [(1, 0),
                (2, 10),
                (3, 20),
                (4, 30),
                (5, 40)]

        window = SlidingWindow(3, mean)

        result, time = window.push(vals[0])
        assert np.isnan(result)
        assert time == 0

        result, time = window.push(vals[1])
        assert np.isnan(result)
        assert time == 10

        result, time = window.push(vals[2])
        assert result == 2
        assert time == 20

        result, time = window.push(vals[3])
        assert result == 3
        assert time == 30

        result, time = window.push(vals[4])
        assert result == 4
        assert time == 40
    
    def test_uses_time_normalized_mean(self):
        vals = [(1, 0),
                (2, 10),
                (3, 20),
                (4, 30),
                (5, 40)]
        window = SlidingWindow(3, normalized_mean, uses_time=True)
        results = [window.push(v) for v in vals]

        back_val, back_time = results[-1]
        assert back_val == 4 / 20
        assert back_time == 40


    def test_window_size_too_small(self):
        with pytest.raises(TypeError):
            window = SlidingWindow(1)
    

class TestSlidingCutoff:

    @pytest.mark.parametrize("cutoff_func", [all_cutoff, median_cutoff])
    def test_window_size_transition(self, cutoff_func):
        ''' Test that saturation is not reached until
        window_size *sliding windows* have been observed.
        '''
        vals = [5, 5, 5, 5]
        cutoff = 2
        window = SlidingCutoff(3, mean, cutoff_func(cutoff))

        results = [window.push(v) for v in vals]
        pprint(results)
        assert results[-1][0] is False
        assert results[-1][1] == 5
        assert results[-1][2] == 3

        reached, smoothed, time = window.push(5)
        assert reached is True
        assert smoothed == 5
        assert time == 4