import pytest

from utils import *
from boink.dbg_tests import (_test_add_single_kmer,
                             _test_add_two_kmers,
                             _test_kmer_degree,
                             _test_kmer_in_degree,
                             _test_kmer_out_degree,
                             _test_kmer_count,
                             _test_add_loop)
from boink.dbg import ExactDBG, kmers


@pytest.fixture
def G(K):
    return ExactDBG(K)

def test_add_single_kmer():
    _test_add_single_kmer()

def test_add_two_kmers():
    _test_add_two_kmers()

def test_kmer_degree():
    _test_kmer_degree()

def test_kmer_in_degree():
    _test_kmer_in_degree()

def test_kmer_out_degree():
    _test_kmer_out_degree()

def test_kmer_count():
    _test_kmer_count()

def test_add_loop(random_sequence, K):
    seq = random_sequence()
    seq += seq[:K]
    _test_add_loop(seq, K)

def test_right_tip(right_tip_structure, G, K):
    (sequence, tip), S = right_tip_structure
    G.add_sequence(sequence)
    
    for kmer in kmers(sequence[1:-1], K):
        assert G.kmer_count(kmer) == 1
        assert G.kmer_degree(kmer) == 2
    assert G.kmer_degree(sequence[:K]) == 1
    assert G.kmer_degree(sequence[-K:]) == 1

    G.add_sequence(tip)
    assert G.kmer_out_degree(sequence[S-K:S]) == 2
    assert G.kmer_in_degree(sequence[S-K:S]) == 1
    assert G.kmer_in_degree(tip[-K:]) == 1
    assert G.kmer_out_degree(tip[S-K:S]) == 2
    for kmer in kmers(tip[1:-2], K):
        assert G.kmer_count(kmer) == 2
        assert G.kmer_degree(kmer) == 2
