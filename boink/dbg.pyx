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

from lemon cimport (SmartDigraph, Node, NodeIt, NodeMap, CrossRefMap, INVALID,
                    countNodes, countArcs, countOutArcs, countInArcs)
from utils cimport _bstring, _ustring


def kmers(sequence, K):
    for i in range(len(sequence)-K+1):
        yield sequence[i:i+K]


cdef class NodeView:

    def __cinit__(self):
        pass

    @staticmethod
    cdef NodeView wrap(Node node, string kmer, uint8_t in_degree,
            uint8_t out_degree, uint16_t count):

        cdef NodeView nv = NodeView()
        nv._this.reset(new Node(node))
        nv.kmer = kmer
        nv.in_degree = in_degree
        nv.out_degree = out_degree
        nv.count = count

        return nv


cdef class ExactDBG:
    '''An exact de Bruijn graph built with LEMON.
    '''

    def __cinit__(self, K, alphabet='DNA_SIMPLE'):
        self.K = K
        self.alphabet = alphabet
        self._alphabet = Alphabets._get(_bstring(alphabet))
        self.graph = new SmartDigraph()
        self.kmers = new CrossRefMap[SmartDigraph,Node,string](deref(self.graph))
        self.counts = new NodeMap[uint16_t](deref(self.graph))

    def nodes(self):
        cdef NodeIt n = NodeIt(deref(self.graph))
        while n != INVALID:
            yield NodeView.wrap(<Node>n, deref(self.kmers)[<Node>n], 
                                self._node_in_degree(<Node>n),
                                self._node_out_degree(<Node>n), 
                                self._node_count(<Node>n))
            preinc(n)

    cdef uint64_t _n_nodes(self):
        return countNodes[SmartDigraph](deref(self.graph))

    @property
    def n_nodes(self):
        return self._n_nodes()

    cdef uint64_t _n_arcs(self):
        return countArcs[SmartDigraph](deref(self.graph))

    @property
    def n_arcs(self):
        return self._n_arcs()

    @property
    def n_edges(self):
        return self._n_arcs()

    cdef uint8_t _node_in_degree(self, Node node):
        return countInArcs[SmartDigraph](deref(self.graph), node)

    cdef uint8_t _kmer_in_degree(self, string kmer):
        return self._node_in_degree(deref(self.kmers).get(kmer))

    def kmer_in_degree(self, kmer):
        cdef string cpkmer = _bstring(kmer)
        return self._kmer_in_degree(cpkmer)

    cdef uint8_t _node_out_degree(self, Node node):
        return countOutArcs[SmartDigraph](deref(self.graph), node)

    cdef uint8_t _kmer_out_degree(self, string kmer):
        return self._node_out_degree(deref(self.kmers).get(kmer))

    def kmer_out_degree(self, kmer):
        cdef string cpkmer = _bstring(kmer)
        return self._kmer_out_degree(cpkmer)

    cdef uint8_t _node_degree(self, Node node):
        return self._node_in_degree(node) + self._node_out_degree(node)

    cdef uint8_t _kmer_degree(self, string kmer):
        return self._node_degree(deref(self.kmers).get(kmer))

    def kmer_degree(self, kmer):
        cdef string cpkmer = _bstring(kmer)
        return self._kmer_degree(cpkmer)

    cdef uint16_t _node_count(self, Node node):
        return deref(self.counts)[node]

    cdef uint16_t _kmer_count(self, string kmer):
        return self._node_count(deref(self.kmers).get(kmer))

    def kmer_count(self, kmer):
        cdef string cpkmer = _bstring(kmer)
        return self._kmer_count(cpkmer)

    cdef uint16_t _inc_node_count(self, Node node):
        cdef uint16_t count = self._node_count(node) + 1
        deref(self.counts).set(node, count)
        return count

    cdef uint16_t _inc_kmer_count(self, string kmer):
        return self._inc_node_count(deref(self.kmers).get(kmer))

    def inc_kmer_count(self, kmer):
        cdef string cpkmer = _bstring(kmer)
        return self._inc_kmer_count(cpkmer)

    def add_sequence(self, sequence):
        cdef string cpsequence = _bstring(sequence)
        return self._add_sequence(cpsequence)
    consume = add_sequence # alias for oxli compat

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
    count = add_kmer # alias for oxli compat

    @cython.boundscheck(False)
    @cython.nonecheck(False)
    @cython.wraparound(False)
    cdef int _add_kmer(self, string kmer) except -1:
        if kmer.length() != self.K:
            raise ValueError('Sequence must have length K={0}; got {1}'.format(self.K, _ustring(kmer)))
        cdef Node node = deref(self.kmers).get(kmer)
        cdef uint16_t count
        #print('k-mer:', _ustring(kmer))
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
        for b in self._alphabet:
            neighbor_kmer[0] = b
            #print('<--', _ustring(neighbor_kmer))
            neighbor_node = deref(self.kmers).get(neighbor_kmer)
            if neighbor_node != INVALID:
                #print('** neighbor:', _ustring(neighbor_kmer))
                deref(self.graph).addArc(neighbor_node, node)

        # copy the root suffix into the neighbor prefix
        # check for right neighbor
        for i in range(1, self.K):
            neighbor_kmer[i-1] = kmer[i]
        for b in self._alphabet:
            neighbor_kmer[self.K-1] = b
            #print('-->', _ustring(neighbor_kmer))
            neighbor_node = deref(self.kmers).get(neighbor_kmer)
            if neighbor_node != INVALID:
                #print('** neighbor:', _ustring(neighbor_kmer))
                deref(self.graph).addArc(node, neighbor_node)

        return True 
        
