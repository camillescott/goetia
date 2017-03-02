# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8, embedsignature=True

import cython

from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc

from libcpp cimport bool
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr
from libcpp.string cimport string
from libc.stdint cimport uint8_t, uint16_t, uint64_t

from khmer._oxli.wrapper cimport _revcomp
from khmer._oxli.parsing cimport Alphabets

from lemon cimport (SmartDigraph, Node, NodeMap, CrossRefMap, INVALID,
                    countNodes, countArcs, countOutArcs, countInArcs)
from utils cimport _bstring, _ustring


cdef class ExactDBG:
    '''An exact de Bruijn graph built with LEMON.
    '''

    def __cinit__(self, K, alphabet='DNA_SIMPLE'):
        self.K = K
        self.alphabet = Alphabets._get(_bstring(alphabet))
        self.graph = new SmartDigraph()
        self.kmers = new CrossRefMap[SmartDigraph,Node,string](deref(self.graph))
        self.counts = new NodeMap[uint16_t](deref(self.graph))

    cdef uint64_t _n_nodes(self):
        return countNodes[SmartDigraph](deref(self.graph))

    cdef uint64_t _n_arcs(self):
        return countArcs[SmartDigraph](deref(self.graph))

    cdef uint8_t _node_in_degree(self, Node node):
        return countInArcs[SmartDigraph](deref(self.graph), node)

    cdef uint8_t _kmer_in_degree(self, string kmer):
        return self._node_in_degree(deref(self.kmers).get(kmer))

    cdef uint8_t _node_out_degree(self, Node node):
        return countOutArcs[SmartDigraph](deref(self.graph), node)

    cdef uint8_t _kmer_out_degree(self, string kmer):
        return self._node_out_degree(deref(self.kmers).get(kmer))

    cdef uint8_t _node_degree(self, Node node):
        return self._node_in_degree(node) + self._node_out_degree(node)

    cdef uint8_t _kmer_degree(self, string kmer):
        return self._node_degree(deref(self.kmers).get(kmer))

    cdef uint16_t _node_count(self, Node node):
        return deref(self.counts)[node]

    cdef uint16_t _kmer_count(self, string kmer):
        return self._node_count(deref(self.kmers).get(kmer))

    cdef uint16_t _inc_node_count(self, Node node):
        cdef uint16_t count = self._node_count(node) + 1
        deref(self.counts).set(node, count)
        return count

    cdef uint16_t _inc_kmer_count(self, string kmer):
        return self._inc_node_count(deref(self.kmers).get(kmer))

    cdef int _add_sequence(self, string sequence) except -1:
        if sequence.length() < self.K:
            raise ValueError('Sequence too short.')
        cdef int i
        cdef int n_new_kmers = 0
        for i in range(0,sequence.length()-self.K+1):
            n_new_kmers += self._add_kmer(sequence.substr(i,self.K))
        return n_new_kmers

    def add_kmer(self, kmer):
        cdef string cpkmer = _bstring(kmer)
        return self._add_kmer(cpkmer)

    @cython.boundscheck(False)
    @cython.nonecheck(False)
    @cython.wraparound(False)
    cdef int _add_kmer(self, string kmer) except -1:
        if kmer.length() != self.K:
            raise ValueError('Sequence must have length K={0}; got {1}'.format(self.K, _ustring(kmer)))
        cdef Node node = deref(self.kmers).get(kmer)
        cdef uint16_t count
        if node != INVALID:
            self._inc_node_count(node)
            return False

        node = deref(self.graph).addNode()
        deref(self.counts).set(node, 1)
        deref(self.kmers).set(node, kmer)
        
        cdef string neighbor_kmer = string(kmer)
        cdef Node neighbor_node
        cdef char  b
        cdef int i

        # copy the root prefix into the neighbor suffix
        # check for left neighbor
        for i in range(0,self.K-1):
            neighbor_kmer[i+1] = kmer[i]
        for b in self.alphabet:
            neighbor_kmer[0] = b
            #print('<--', _ustring(neighbor_kmer))
            neighbor_node = deref(self.kmers).get(neighbor_kmer)
            if neighbor_node != INVALID:
                deref(self.graph).addArc(neighbor_node, node)

        # copy the root suffix into the neighbor prefix
        # check for right neighbor
        for i in range(1, self.K):
            neighbor_kmer[i-1] = kmer[i]
        for b in self.alphabet:
            neighbor_kmer[self.K-1] = b
            #print('-->', _ustring(neighbor_kmer))
            neighbor_node = deref(self.kmers).get(neighbor_kmer)
            if neighbor_node != INVALID:
                deref(self.graph).addArc(node, neighbor_node)

        return True 
        
