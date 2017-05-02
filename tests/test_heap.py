import pytest

import random
from boink.heap import Weighted, MinHeap

@pytest.mark.parametrize('capacity', [10, 100, 1000])
def test_topN_stream(capacity):
    stream = range(2 * capacity)
    heap = MinHeap(capacity)

    for val in stream:
        heap.insert(Weighted(val, val))

    assert sorted(heap.values()) == list(range(capacity, 2*capacity))
    assert heap.n_items == capacity

@pytest.mark.parametrize('capacity', [10])
def test_topN_random(capacity):
    sample_size = capacity * 10
    stream = [random.randint(0, sample_size) for _ in range(sample_size)]
    heap = MinHeap(capacity)

    for val in stream:
        heap.insert(Weighted(val, val))

    assert sorted(heap.values()) == sorted(stream)[-capacity:]
    assert heap.n_items == capacity
