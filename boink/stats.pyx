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

from khmer._oxli.wrapper cimport (CpHashgraph, get_hashgraph_ptr, 
                                  _hash_forward, _hash_murmur, HashIntoType, 
                                  _revhash, BoundedCounterType, _hash,
                                  CpStreamingPartitioner)
from khmer._oxli.partitioning cimport StreamingPartitioner
from khmer._oxli.assembly cimport CompactingAssembler
from khmer._oxli.hashing cimport Kmer
from khmer._oxli.hashing import Kmer
from khmer._oxli.parsing cimport Sequence
from khmer._oxli.parsing import Sequence
from khmer._oxli.traversal cimport Traverser

from .utils cimport _bstring, _ustring


cdef class GraphFunction:

    def __cinit__(self, object graph=None, *args, **kwargs):
        if graph is not None:
            self.graph = get_hashgraph_ptr(graph)
            self.K = deref(self.graph).ksize()

    cdef void _set_graph(self, CpHashgraph * ptr):
        self.graph = ptr
        self.K = deref(self.graph).ksize()

    def set_graph(self, object graph):
        self.graph = get_hashgraph_ptr(graph)
        self.K = deref(self.graph).ksize()

    cpdef int consume(self, Sequence sequence):
        return deref(self.graph).consume_string(sequence._obj.sequence)

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

    def __cinit__(self, object graph=None, StreamingPartitioner partitioner=None, 
                  *args, **kwargs):
        self.partitioner = partitioner._this
        self.graph = deref(self.partitioner).graph
        self.K = deref(self.graph).ksize()

    cdef float _evaluate_tags(self, string sequence, set[HashIntoType]& tags):
        pass

    cpdef float evaluate_tags(self, Sequence sequence, tuple tags):
        pass


cdef class PartitionCoverage(PartitionFunction):

    def __cinit__(self, coverage_cutoff=20, object graph=None, 
                  StreamingPartitioner partitioner=None, *args, **kwargs):
        self.coverage_cutoff = coverage_cutoff

    cdef float _evaluate_tags(self, string sequence, set[HashIntoType]& tags):
        cdef uint64_t n_tags = len(tags)
        cdef float acc = 0
        cdef uint64_t tag
        for tag in tags:
            acc += <float>deref(self.graph).get_count(tag)
        if acc / <float>n_tags > self.coverage_cutoff:
            return 1.0 # dummy fitness vals
        else:
            return 0.0


cdef class KmerCountFunction(GraphFunction):

    cpdef float evaluate_kmer(self, Kmer kmer):
        return <int>deref(self.graph).get_count(deref(kmer._this.get()).kmer_u)


cdef class KmerDegreeFunction(GraphFunction):

    def __cinit__(self, object graph=None, *args, **kwargs):
        self.traverser = Traverser(graph)

    def set_graph(self, object graph):
        self.graph = get_hashgraph_ptr(graph)
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

    def __cinit__(self, object graph=None, int bandwidth=12, *args, **kwargs):
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
                count = deref(self.graph).get_count(_hash(_bstring(kmer), self.K))
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
