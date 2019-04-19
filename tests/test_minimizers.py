import pytest

from boink.minimizers import InteriorMinimizer

class TestInteriorMinimizer(object):

    def test_basic(self):
        sequence = [3, 1, 2, 7, 8, 4, 6, 8]
        window_size = 3
        expected = [1, 1, 2, 4, 4, 4]

        M = InteriorMinimizer(window_size)
        minimizers = M(sequence)

        assert minimizers == expected

    def test_window_too_large(self):
        sequence = [3, 1, 2, 3, 1]
        M = InteriorMinimizer(6)

        assert M(sequence) == []

    def test_window_equals_seqlen(self):
        sequence = [3, 2, 1, 2, 3]
        M = InteriorMinimizer(5)

        assert M(sequence) == [1]
