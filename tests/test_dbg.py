import pytest

from utils import *
from boink.dbg_tests import (_test_add_kplus1mer,
                             _test_kmer_degree,
                             _test_kmer_in_degree,
                             _test_kmer_out_degree,
                             _test_get_count,
                             _test_add_loop)
from boink.dbg import ExactDBG, kmers


@pytest.fixture
def G(K):
    return ExactDBG(K)

def test_add_kplus1mer():
    _test_add_kplus1mer()


def test_kmer_degree():
    _test_kmer_degree()

def test_kmer_in_degree():
    _test_kmer_in_degree()

def test_kmer_out_degree():
    _test_kmer_out_degree()

def test_get_count():
    _test_get_count()

def test_add_loop(random_sequence, K):
    seq = random_sequence()
    seq += seq[:K]
    _test_add_loop(seq, K)

def test_right_tip(right_tip_subgraph, G, K):
    (sequence, tip), S = right_tip_subgraph
    G.add_sequence(sequence)
    
    for kmer in kmers(sequence[1:-1], K):
        assert G.kmer_degree(kmer) == 2
    for kmer in kmers(sequence[1:-1], K+1):
        assert G.sequence_count(kmer) == 1
    assert G.kmer_degree(sequence[:K]) == 1
    assert G.kmer_degree(sequence[-K:]) == 1

    G.add_sequence(tip)
    assert G.kmer_out_degree(sequence[S-K:S]) == 2
    assert G.kmer_in_degree(sequence[S-K:S]) == 1
    assert G.kmer_in_degree(tip[-K:]) == 1
    assert G.kmer_out_degree(tip[S-K:S]) == 2
    for kmer in kmers(tip[1:-2], K):
        assert G.kmer_degree(kmer) == 2
    for kmer in kmers(tip[1:-2], K+1):
        assert G.sequence_count(kmer) == 2

def test_left_tip(left_tip_subgraph, G, K):
    (sequence, tip), S = left_tip_subgraph
    G.add_sequence(sequence)

    for kmer in kmers(sequence[1:-1], K):
        assert G.kmer_degree(kmer) == 2
    for kmer in kmers(sequence[1:-1], K+1):
        assert G.sequence_count(kmer) == 1
    assert G.kmer_degree(sequence[:K]) == 1
    assert G.kmer_degree(sequence[-K:]) == 1

    G.add_sequence(tip)
    assert G.kmer_in_degree(sequence[S+1:S+1+K]) == 2
    assert G.kmer_in_degree(tip[1:K+1]) == 2
    assert G.kmer_out_degree(tip[0:K]) == 1
    assert G.kmer_in_degree(tip[0:K]) == 0

    for kmer in kmers(tip[1:], K+1):
        assert G.sequence_count(kmer) == 2

def test_snp_bubble(snp_bubble_subgraph, G, K):
    (wildtype, mutant), S = snp_bubble_subgraph
    G.add_sequence(wildtype)
    for kmer in kmers(wildtype[1:-1], K):
        #assert G.kmer_count(kmer) == 1
        assert G.kmer_degree(kmer) == 2
    assert G.kmer_degree(wildtype[:K]) == 1
    assert G.kmer_degree(wildtype[-K:]) == 1

    G.add_sequence(mutant)
    assert G.kmer_degree(wildtype[:K]) == 1
    assert G.kmer_degree(wildtype[-K:]) == 1
    #assert G.kmer_count(wildtype[:K]) == 2
    #assert G.kmer_count(wildtype[-K:]) == 2
    
    # Check the left HDN
    assert G.kmer_out_degree(wildtype[S-K:S]) == 2
    assert G.kmer_in_degree(wildtype[S-K:S]) == 1

    # Check the wildtype bubble path
    for kmer in kmers(wildtype[S-K+1:S+K], K):
        assert G.kmer_degree(kmer) == 2
        #assert G.kmer_count(kmer) == 1

    # Check the mutant bubble path
    for kmer in kmers(mutant[S-K+1:S+K], K):
        assert G.kmer_degree(kmer) == 2
        #assert G.kmer_count(kmer) == 1

    # Check the right HDN
    assert G.kmer_out_degree(wildtype[S+1:S+K+1]) == 1
    assert G.kmer_in_degree(wildtype[S+1:S+K+1]) == 2


