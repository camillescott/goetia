# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8

from libcpp cimport bool
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint16_t, uint64_t

cimport numpy as np

from lemon cimport SmartDigraph, NodeMap, Node, INVALID, CrossRefMap
from khmer._oxli.wrapper cimport _revcomp, CpCounttable
from khmer._oxli.parsing cimport Alphabets



cdef class NodeView:

    cdef shared_ptr[Node] _this
    cdef public string kmer
    cdef public uint8_t in_degree
    cdef public uint8_t out_degree

    @staticmethod
    cdef NodeView wrap(Node node, string kmer, uint8_t in_degree,
                       uint8_t out_degree)


cdef class ExactDBG:

    cdef SmartDigraph * graph
    cdef CrossRefMap[SmartDigraph, Node, string] * kmers
    #cdef NodeMap[uint16_t] * counts
    cdef CpCounttable * counts

    cdef readonly int K
    cdef readonly object alphabet
    cdef string _alphabet

    cdef int _add_sequence(self, string) except -1
    cdef int _add_kplus1mer(self, string) except -1
    cdef int _consume_fastx(self, unicode filename) except -1

    cdef uint64_t _n_nodes(self)
    cdef uint64_t _n_arcs(self)

    cdef uint8_t _node_in_degree(self, Node node)
    cdef uint8_t _kmer_in_degree(self, string kmer)

    cdef uint8_t _node_out_degree(self, Node node)
    cdef uint8_t _kmer_out_degree(self, string kmer)

    cdef uint8_t _node_degree(self, Node node)
    cdef uint8_t _kmer_degree(self, string kmer)

    cdef np.ndarray[np.uint16_t] _kmer_counts(self, string kmer)

    cdef bool _inc_count(self, string sequence)
    cdef uint16_t _get_count(self, string sequence)


cdef class CompactingDBG(ExactDBG):

    #cdef void _gather_linear_path(self, vector[string]& path)
    #cdef string _extract_contig(self, vector[string]& path)
    pass
