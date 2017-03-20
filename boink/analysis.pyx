# cython: c_string_type=unicode, c_string_encoding=utf8, embedsignature=True

from cython.operator cimport dereference as deref

from khmer import Nodegraph
from khmer._oxli.parsing cimport SanitizedFastxParser, Sequence
from khmer._oxli.hashing cimport Kmer, CpKmer
from khmer._oxli.wrapper cimport CpKmerIterator, CpNodegraph, CpCountgraph
from khmer._oxli.wrapper cimport get_hashgraph_ptr
from .heap cimport MinHeap, Weighted
from .stats cimport GraphFunction, KmerCountFunction, KmerDegreeFunction


def find_top_N_kmers(list filenames, int N, GraphFunction func,
                     object graph):

    cdef object seen = Nodegraph(graph.ksize(), max(graph.hashsizes()), 
                                 len(graph.hashsizes()))
    cdef MinHeap minheap = MinHeap(N)
    cdef unicode filename
    for filename in filenames:
        print('Processing {0}'.format(filename))
        _find_top_N_kmers(filename, N, func, <CpCountgraph *>get_hashgraph_ptr(graph),
                <CpNodegraph *>get_hashgraph_ptr(seen), minheap)
    cdef dict results = dict(minheap.iteritems())

    return results


cdef _find_top_N_kmers(unicode filename, int N, GraphFunction func,
                       CpCountgraph * counts, CpNodegraph * seen, 
                       MinHeap minheap):

    cdef SanitizedFastxParser parser = SanitizedFastxParser(filename)
    cdef Sequence sequence
    cdef int K = deref(counts).ksize()
    cdef CpKmer cpkmer
    cdef Kmer kmer
    cdef CpKmerIterator * it
    cdef int n_kmers = 0

    cdef float weight

    for n_reads, sequence in enumerate(parser):
        if n_reads % 50000 == 0:
            print('Checked {0} reads, {1} k-mers...'.format(n_reads, n_kmers))
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

def find_N_most_abundant_kmers(filenames, N, graph):
    return find_top_N_kmers(filenames, N, KmerCountFunction(graph), graph)
