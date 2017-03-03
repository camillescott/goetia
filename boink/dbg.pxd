# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8

from libcpp cimport bool
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr
from libcpp.string cimport string
from libc.stdint cimport uint8_t, uint16_t, uint64_t

from lemon cimport SmartDigraph, NodeMap, Node, INVALID, CrossRefMap
from khmer._oxli.wrapper cimport _revcomp
from khmer._oxli.parsing cimport Alphabets


cdef class NodeView:

    cdef shared_ptr[Node] _this
    cdef public string kmer
    cdef public uint8_t in_degree
    cdef public uint8_t out_degree
    cdef public uint16_t count

    @staticmethod
    cdef NodeView wrap(Node node, string kmer, uint8_t in_degree,
                       uint8_t out_degree, uint16_t count)


cdef class ExactDBG:

    cdef SmartDigraph * graph
    cdef CrossRefMap[SmartDigraph, Node, string] * kmers
    cdef NodeMap[uint16_t] * counts

    cdef readonly int K
    cdef readonly object alphabet
    cdef string _alphabet

    cdef int _add_sequence(self, string) except -1
    cdef int _add_kmer(self, string) except -1

    cdef uint64_t _n_nodes(self)
    cdef uint64_t _n_arcs(self)

    cdef uint8_t _node_in_degree(self, Node node)
    cdef uint8_t _kmer_in_degree(self, string kmer)
    cdef uint8_t _node_out_degree(self, Node node)
    cdef uint8_t _kmer_out_degree(self, string kmer)
    cdef uint8_t _node_degree(self, Node node)
    cdef uint8_t _kmer_degree(self, string kmer)
    cdef uint16_t _node_count(self, Node node)
    cdef uint16_t _kmer_count(self, string kmer)
    cdef uint16_t _inc_node_count(self, Node node)
    cdef uint16_t _inc_kmer_count(self, string kmer)

