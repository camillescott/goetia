import pytest

from boink.dbg_tests import (_test_add_single_kmer,
                             _test_add_two_kmers,
                             _test_kmer_degree,
                             _test_kmer_in_degree,
                             _test_kmer_out_degree)

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
