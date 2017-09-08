from khmer._oxli.oxli_types cimport *
from khmer._oxli.graphs cimport CpHashgraph, Hashgraph
from khmer._oxli.hashing cimport CpKmer, Kmer
from khmer._oxli.parsing cimport Sequence

from libc.stdint cimport uint32_t, uint64_t
from libcpp cimport bool

from stats cimport GraphFunction

cdef class ProbabilisticGraph:

    cdef Hashgraph * graph 
    cdef WordLength K
    cdef GraphFunction func

    cdef bool _insert_kmer(self, Kmer kmer)
    cdef bool _insert_sequence(self, Sequence sequence)

