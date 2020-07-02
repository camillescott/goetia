# goetia/tests/test_saturation.py
# Copyright (C) 2020 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from statistics import mean

import numpy as np
import pytest

from goetia.saturation import SlidingWindow, SlidingCutoff, normalized_mean


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
    
    