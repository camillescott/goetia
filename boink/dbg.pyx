# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8, embedsignature=True
from __future__ import print_function

import cython
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc

cimport numpy as np
import sys

from libcpp cimport bool
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr
from libcpp.string cimport string
from libc.stdint cimport uint8_t, uint16_t, uint64_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.hashing cimport _revcomp, _hash_murmur
from khmer._oxli.parsing cimport Alphabets, Sequence, SanitizedFastxParser
from khmer._oxli.utils cimport _get_n_primes_near_x

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
    cdef NodeView wrap(Node node, string kmer, 
                       uint8_t in_degree, uint8_t out_degree):

        cdef NodeView nv = NodeView()
        nv._this.reset(new Node(node))
        nv.kmer = kmer
        nv.in_degree = in_degree
        nv.out_degree = out_degree

        return nv


cdef class ExactDBG:
    '''An exact de Bruijn graph built with LEMON.
    '''

    def __cinit__(self, K, alphabet='DNA_SIMPLE', n_tables=4, table_size=1e6):
        self.K = K
        self.alphabet = alphabet
        self._alphabet = Alphabets._get(_bstring(alphabet))
        self.graph = new SmartDigraph()
        self.kmers = new CrossRefMap[SmartDigraph,Node,string](deref(self.graph))
        #self.counts = new NodeMap[uint16_t](deref(self.graph))
        self.counts = new CpCounttable(K+1, 
                          _get_n_primes_near_x(n_tables, table_size))

    def nodes(self):
        cdef NodeIt n = NodeIt(deref(self.graph))
        while n != INVALID:
            yield NodeView.wrap(<Node>n, deref(self.kmers)[<Node>n], 
                                self._node_in_degree(<Node>n),
                                self._node_out_degree(<Node>n))
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


    def sequence_count(self, sequence):
        cdef string cseq = _bstring(sequence)
        return self._get_count(cseq)

    cdef uint16_t _get_count(self, string sequence):
        cdef HashIntoType h = _hash_murmur(sequence, sequence.length())
        return deref(self.counts).get_count(h)

    cdef bool _inc_count(self, string sequence):
        cdef HashIntoType h = _hash_murmur(sequence, sequence.length())
        return deref(self.counts).add(h)


    def add_sequence(self, sequence):
        cdef string cpsequence = _bstring(sequence)
        return self._add_sequence(cpsequence)
    consume = add_sequence # alias for oxli compat

    cdef int _add_sequence(self, string sequence) except -1:
        if sequence.length() < self.K+1:
            raise ValueError('Sequence too short.')
        cdef int i
        cdef int n_new_kmers = 0
        for i in range(0,sequence.length()-self.K):
            n_new_kmers += self._add_kplus1mer(sequence.substr(i,self.K+1))
        return n_new_kmers

    cdef np.ndarray[np.uint16_t] _kmer_counts(self, string sequence):
        if sequence.length() < self.K+1:
            raise ValueError('Sequence too short.')
        cdef int i
        cdef np.ndarray[np.uint16_t] counts = np.zeros([sequence.length()-self.K],
                                                      dtype=np.uint16)
        for i in range(0,sequence.length()-self.K):
            counts[i] = self._kmer_count(sequence.substr(i,self.K+1))

        return counts

    def add_kplus1mer(self, kp1mer):
        if len(kp1mer) != self.K+1:
            raise ValueError('Sequence must have length K={0}; got'\
                    '{1}'.format(self.K+1, kp1mer))

        cdef string cpkp1mer = _bstring(kp1mer)
        return self._add_kplus1mer(cpkp1mer)
    count = add_kplus1mer # alias for oxli compat

    @cython.boundscheck(False)
    @cython.nonecheck(False)
    @cython.wraparound(False)
    cdef int _add_kplus1mer(self, string kp1mer) except -1:

        if self._get_count(kp1mer):
            self._inc_count(kp1mer)
            return 0
            
        cdef int n_new = 0
        cdef Node left = deref(self.kmers).get(kp1mer.substr(0, self.K))
        cdef Node right = deref(self.kmers).get(kp1mer.substr(1, self.K))

        if left == INVALID:
            left = deref(self.graph).addNode()
            deref(self.kmers).set(left, kp1mer.substr(0, self.K))
            n_new += 1

        if right == INVALID:
            right = deref(self.graph).addNode()
            deref(self.kmers).set(right, kp1mer.substr(1, self.K))
            n_new += 1

        deref(self.graph).addArc(left, right)
        self._inc_count(kp1mer)

        return n_new

    cdef int _consume_fastx(self, unicode filename) except -1:
        cdef SanitizedFastxParser parser = SanitizedFastxParser(filename)
        cdef Sequence sequence
        cdef int n_reads
        cdef int n_kmers = 0

        for n_reads, sequence in enumerate(parser):
            if n_reads % 100000 == 0:
                print('Consumed {0} reads...'.format(n_reads), file=sys.stderr)
            n_kmers += self._add_sequence(sequence._obj.sequence)

        return n_kmers

    def consume_fastx(self, filename):
        return self._consume_fastx(filename)
