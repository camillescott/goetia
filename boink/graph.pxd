from khmer._oxli._oxli cimport CpHashgraph, CpKmer, WordLength
from khmer._oxli.hashing cimport Kmer

from libcpp.vector cimport vector
from libc.stdint cimport uint32_t, uint64_t
from libcpp.memory import unique_ptr, shared_ptr

cdef class ProbabilisticGraph:

    cdef shared_ptr[CpHashgraph] _graph
    cdef likelihood(self, CpKmer)

    cdef WordLength K
    cdef vector[uint64_t] table_sizes
    cdef uint32_t n_tables
