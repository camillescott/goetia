# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8
from __future__ import print_function
import sys
from libcpp.memory cimport shared_ptr, weak_ptr
from libcpp.string cimport string
from libcpp.set cimport set
cimport cython
import cython
from cython.operator cimport dereference as deref
from libc.stdint cimport uint8_t, uint16_t, uint64_t

from khmer._oxli.hashing cimport (_hash_forward, _hash_murmur,
                                  _revhash, _hash)
from khmer._oxli.graphs cimport CpHashgraph, Hashgraph
from khmer._oxli.partitioning cimport (StreamingPartitioner,
                                       CpStreamingPartitioner)
from khmer._oxli.assembly cimport CompactingAssembler
from khmer._oxli.hashing cimport Kmer
from khmer._oxli.sequence cimport Sequence
from khmer._oxli.traversal cimport Traverser

from .utils cimport _bstring, _ustring


cdef class GraphFunction:

    def __cinit__(self, Hashgraph graph=None, *args, **kwargs):
        if graph is not None:
            self.graph = graph
            self.K = self.graph.ksize()

    def set_graph(self, Hashgraph graph not None):
        self.graph = graph
        self.K = self.graph.ksize()

    cpdef int consume(self, Sequence sequence):
        return deref(self.graph._hg_this).consume_string(sequence._obj.sequence)

    cpdef float evaluate_kmer(self, Kmer kmer):
        raise NotImplementedError()

    cpdef float evaluate_sequence(self, Sequence sequence):
        raise NotImplementedError()

    def save(self, filename):
        raise NotImplementedError()

    @staticmethod
    def load(filename):
        raise NotImplementedError()


cdef class PartitionFunction(GraphFunction):

    def __cinit__(self, Hashgraph graph=None, StreamingPartitioner partitioner=None, 
                  *args, **kwargs):
        self.partitioner = partitioner._this
        if partitioner is not None:
            self.graph = partitioner.graph
        self.K = self.graph.ksize()

    cdef float _evaluate_tags(self, string sequence, vector[HashIntoType]& tags):
        pass

    cpdef float evaluate_tags(self, Sequence sequence, tuple tags):
        pass


cdef class PartitionCoverage(PartitionFunction):

    def __cinit__(self, coverage_cutoff=20, object graph=None, 
                  StreamingPartitioner partitioner=None, *args, **kwargs):
        self.coverage_cutoff = coverage_cutoff

    cdef float _evaluate_tags(self, string sequence, vector[HashIntoType]& tags):
        cdef uint64_t n_tags = len(tags)
        cdef float acc = 0
        cdef uint64_t tag
        for tag in tags:
            acc += <float>self.graph.get_count(tag)
        cdef float val = (acc / <float>n_tags)
        return val / <float>self.coverage_cutoff


cdef class PartitionCoverageSlice(PartitionFunction):

    def __cinit__(self, floor=10, ceiling=15, object graph=None, 
                  StreamingPartitioner partitioner=None, *args, **kwargs):
        self.coverage_floor = floor
        self.coverage_ceiling = ceiling

    cdef float _evaluate_tags(self, string sequence, vector[HashIntoType]& tags):
        cdef uint64_t n_tags = len(tags)
        cdef float acc = 0
        cdef uint64_t tag
        for tag in tags:
            acc += <float>self.graph.get_count(tag)
        cdef float val = (acc / <float>n_tags)
        if self.coverage_floor < val < self.coverage_ceiling:
            return 1.0
        else:
            return 0.0


cdef class KmerCountFunction(GraphFunction):

    cpdef float evaluate_kmer(self, Kmer kmer):
        return <int>self.graph.get_count(kmer)


cdef class KmerDegreeFunction(GraphFunction):

    def __cinit__(self, Hashgraph graph=None, *args, **kwargs):
        self.traverser = Traverser(graph)

    def set_graph(self, Hashgraph graph):
        self.graph = graph
        self.traverser = Traverser(graph)

    cpdef float evaluate_kmer(self, Kmer kmer):
        return <int>self.traverser.degree(kmer)


cdef class KmerDegreeBiasFunction(KmerDegreeFunction):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef float evaluate_kmer(self, Kmer kmer):
        cdef int degree = self.traverser.degree(kmer)
        self.degree_counts[degree] += 1
        cdef int i
        cdef string kmer_s = _bstring(kmer.kmer)
        cdef HashIntoType h
        for i in range(0, self.K-6):
            h = _hash_forward(kmer_s.substr(i, 7).c_str(), 7)
            self.heptamer_counts[h][degree] += 1
        return degree

    def save(self, fp):
        cdef int i
        cdef str kmer
        try:
            output = open(fp, 'w')
        except TypeError:
            output = fp
        finally:
            for i in range(N_HEPTAMERS):
                row = [str(d) for d in self.heptamer_counts[i][:MAX_DBG_DEGREE]]
                kmer = _revhash(<HashIntoType>i, 7)
                fp.write(kmer)
                fp.write(', ')
                fp.write(', '.join(row))
                fp.write('\n')


cdef class SamplePathsFunction(GraphFunction):

    def __cinit__(self, Hashgraph graph=None, int bandwidth=12, *args, **kwargs):
        self.assembler = CompactingAssembler(graph)
        self.bandwidth = <uint64_t>bandwidth
        self.n_batches = 2 ** self.bandwidth
        self.results = []
        self.t = 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef float evaluate_sequence(self, Sequence sequence):
        cdef unicode kmer
        cdef HashIntoType h
        cdef int length
        cdef BoundedCounterType count
        for kmer in sequence.kmers(self.K):
            #kmer = _ustring(kmer)
            h = _hash_murmur(_bstring(kmer), self.K)
            if (h & (self.n_batches-1)) == 0:
                length = len(self.assembler.assemble(kmer))
                count = self.graph.get_count(kmer)
                #print(h, kmer, file=sys.stderr)
                self.results.append((self.t, h, kmer, length, count))
                if len(self.results) % 1000 == 0:
                    print('{0} k-mers...'.format(len(self.results)), file=sys.stderr)
        self.t += 1

    def save(self, fp):
        cdef unicode kmer
        cdef int length
        cdef BoundedCounterType count
        try:
            output = open(fp, 'w')
        except TypeError:
            output = fp
        finally:
            output.write(','.join(['t', 'hash', 'kmer', 'path_len', 'kmer_count']))
            output.write('\n')
            for row in self.results:
                output.write(','.join((str(item) for item in row)))
                output.write('\n')


FUNCTIONS = [name for name in dir() if name.endswith('Function') \
             and globals()[name] != GraphFunction]
