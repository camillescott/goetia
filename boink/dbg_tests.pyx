# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8

from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr
from libcpp.string cimport string

from lemon cimport SmartDigraph, Node, NodeIt, Arc, INVALID
from dbg import ExactDBG
from dbg cimport ExactDBG
from utils cimport _bstring


cdef ExactDBG get_test_graph(string s, int K):
    cdef ExactDBG g = ExactDBG(K)
    g._add_sequence(s)
    return g


cpdef _test_add_single_kmer():
    cdef string s = 'AAAA'
    cdef ExactDBG g = get_test_graph(s, 4)
    assert g._n_nodes() == 1
    assert g._kmer_count(s) == 1


cpdef _test_add_two_kmers():
    cdef string s = 'ACCCG'
    cdef ExactDBG g = get_test_graph(s, 4)
    assert g._n_nodes() == 2
    assert g._n_arcs() == 1
    assert g._kmer_count(s.substr(1)) == 1
    assert g._kmer_count(s.substr(0,4)) == 1


cpdef _test_kmer_out_degree():
    cdef string s = 'ACCGT'
    cdef ExactDBG g = get_test_graph(s, 4)
    assert g._kmer_out_degree(s.substr(0,4)) == 1
    assert g._kmer_out_degree(s.substr(1,4)) == 0

cpdef _test_kmer_in_degree():
    cdef string s = 'ACCGT'
    cdef ExactDBG g = get_test_graph(s, 4)
    assert g._kmer_in_degree(s.substr(0,4)) == 0
    assert g._kmer_in_degree(s.substr(1,4)) == 1

cpdef _test_kmer_degree():
    cdef string s = 'ACCGT'
    cdef ExactDBG g = get_test_graph(s, 4)
    assert g._kmer_degree(s.substr(0,4)) == 1
    assert g._kmer_degree(s.substr(1,4)) == 1

cpdef _test_kmer_count():
    cdef string s = 'AAAA'
    cdef ExactDBG g = get_test_graph(s, 4)

    for i in range(1, 10):
        assert g._kmer_count(s) == i
        g._add_kmer(s)

cpdef _test_add_loop(seq, K):
    cdef string s = _bstring(seq)
    cdef ExactDBG g = get_test_graph(s, K)
    #print([n.kmer for n in g.nodes()])

    for i in range(0,s.length()-K+1):
        d = g._kmer_degree(s.substr(i,K))
        #print(s.substr(i,K))
        assert d == 2, (d, g._n_nodes())
