# cython: c_string_type=unicode, c_string_encoding=utf8, embedsignature=True, language_level=3

from cython.operator cimport dereference as deref
from libcpp cimport bool
from libcpp.memory cimport shared_ptr

import sys

from khmer._oxli.graphs cimport Nodegraph, Hashgraph, Countgraph
from khmer._oxli.parsing cimport SanitizedFastxParser, Sequence
from khmer._oxli.hashing cimport Kmer, CpKmer, CpKmerIterator
from khmer._oxli.graphs cimport (Nodegraph, CpNodegraph, 
                                 CpCountgraph)
from .heap cimport MinHeap, Weighted
from .stats cimport (GraphFunction, KmerCountFunction, KmerDegreeFunction,
                     KmerDegreeBiasFunction)


def find_top_N_kmers(list filenames, int N, GraphFunction func,
                     Countgraph graph):

    cdef Nodegraph seen = Nodegraph(graph.ksize(), max(graph.hashsizes()), 
                                 len(graph.hashsizes()))
    cdef MinHeap minheap = MinHeap(N)
    for filename in filenames:
        print('Processing {0}'.format(filename), file=sys.stderr)
        _find_top_N_kmers(filename, N, func,
                          seen._ng_this, minheap)
    cdef dict results = dict(minheap.iteritems())

    return results


cdef _find_top_N_kmers(unicode filename, int N, GraphFunction func,
                       shared_ptr[CpNodegraph] seen, 
                       MinHeap minheap):

    cdef SanitizedFastxParser parser = SanitizedFastxParser(filename)
    cdef Sequence sequence
    cdef int K = deref(seen).ksize()
    cdef CpKmer cpkmer
    cdef Kmer kmer
    cdef CpKmerIterator * it
    cdef int n_kmers = 0

    cdef float weight

    for n_reads, sequence in enumerate(parser):
        if n_reads % 50000 == 0:
            print('Checked {0} reads, {1} k-mers...'.format(n_reads, n_kmers),
                  file=sys.stderr)
        it = new CpKmerIterator(sequence._obj.sequence.c_str(), K)
        while not it.done():
            cpkmer = it.next()
            if deref(seen).get_count(cpkmer.kmer_u):
                continue
            kmer = Kmer.wrap(new CpKmer(cpkmer), K)
            weight = func.evaluate_kmer(kmer)
            minheap.insert(Weighted(kmer, weight))
            deref(seen).add(cpkmer.kmer_u)
            n_kmers += 1
        del it


def kmer_apply(list filenames, GraphFunction func):
    for filename in filenames:
        print('Processing {0}'.format(filename), file=sys.stderr)
        _kmer_apply(filename, func)
    return func


cdef _kmer_apply(unicode filename, GraphFunction func):

    cdef SanitizedFastxParser parser = SanitizedFastxParser(filename)
    cdef Sequence sequence
    cdef int K = func.K
    cdef CpKmer cpkmer
    cdef Kmer kmer
    cdef CpKmerIterator * it
    cdef int n_kmers = 0
    cdef int n_reads

    for n_reads, sequence in enumerate(parser):
        if n_reads % 50000 == 0:
            print('Checked {0} reads, {1} k-mers...'.format(n_reads, n_kmers),
                  file=sys.stderr)
        it = new CpKmerIterator(sequence._obj.sequence.c_str(), K)
        while not it.done():
            cpkmer = it.next()
            kmer = Kmer.wrap(new CpKmer(cpkmer), K)
            func.evaluate_kmer(kmer)
            n_kmers += 1
        del it


def sequence_apply(list filenames, GraphFunction func, consume=False):
    for filename in filenames:
        print('Processing {0}'.format(filename), file=sys.stderr)
        _sequence_apply(filename, func, consume)
    return func


cdef _sequence_apply(unicode filename, GraphFunction func, bool consume):
    cdef SanitizedFastxParser parser = SanitizedFastxParser(filename)
    cdef Sequence sequence
    cdef int K = func.K
    cdef int n_reads

    for n_reads, sequence in enumerate(parser):
        if n_reads % 50000 == 0:
            print('Checked {0} reads...'.format(n_reads), file=sys.stderr)
        if consume:
            func.consume(sequence)
        func.evaluate_sequence(sequence)


def find_N_most_abundant_kmers(filenames, N, Countgraph graph):
    return find_top_N_kmers(filenames, N, KmerCountFunction(graph), graph)


def apply_kmer_degree_bias(filenames, Countgraph graph):
    cdef KmerDegreeBiasFunction func = KmerDegreeBiasFunction(graph)
    kmer_apply(filenames, func)
    return func
