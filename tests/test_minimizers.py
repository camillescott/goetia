import pytest

from boink import libboink


class TestInteriorMinimizer(object):

    def test_basic(self):
        sequence = [3, 1, 2, 7, 8, 4, 6, 8]
        window_size = 3
        expected = [1, 1, 2, 4, 4, 4]

        M = libboink.InteriorMinimizer['uint64_t'](window_size)
        for val in sequence:
            M.update(val)

        assert list(M.get_minimizer_values()) == expected

    def test_window_too_large(self):
        sequence = [3, 1, 2, 3, 1]
        M = libboink.InteriorMinimizer['uint64_t'](6)
        for val in sequence:
            M.update(val)
        assert list(M.get_minimizer_values()) == []

    def test_window_equals_seqlen(self):
        sequence = [3, 2, 1, 2, 3]
        M = libboink.InteriorMinimizer['uint64_t'](5)
        for val in sequence:
            M.update(val)
        assert list(M.get_minimizer_values()) ==  [1]
