from khmer._oxli.wrapper cimport CpHashgraph, CpKmer, WordLength
from khmer._oxli.hashing cimport Kmer
from khmer._oxli.parsing cimport Sequence

from libcpp.vector cimport vector
from libc.stdint cimport uint32_t, uint64_t
from libcpp cimport bool

from stats cimport PKmerFunction, PSequenceFunction

cdef class ProbabilisticGraph:

    cdef CpHashgraph * _graph 
    cdef WordLength K
    cdef PKmerFunction kmer_func
    cdef PSequenceFunction sequence_func

    cdef bool _insert_kmer(self, Kmer kmer)
    cdef bool _insert_sequence(self, Sequence sequence)

