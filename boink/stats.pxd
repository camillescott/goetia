from libcpp.memory cimport shared_ptr, weak_ptr

from khmer._oxli._oxli cimport CpHashgraph
from khmer._oxli.hashing cimport Kmer
from khmer._oxli.parsing cimport Sequence


cdef class PKmerFunction:
    cdef float evaluate(self, Kmer kmer, CpHashgraph * graph)


cdef class PSequenceFunction:
    cdef float evaluate(self, Sequence sequence, CpHashgraph * graph)

