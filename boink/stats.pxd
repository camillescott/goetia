from libcpp.memory cimport shared_ptr, weak_ptr

from khmer._oxli.wrapper cimport CpHashgraph
from khmer._oxli.hashing cimport Kmer
from khmer._oxli.parsing cimport Sequence


cdef class PFunction:
    # have PFunction store Hashgraph ptr
    cdef CpHashgraph * graph
    cdef void _set_graph(self, CpHashgraph *)


cdef class PKmerFunction(PFunction):
    cpdef float evaluate(self, Kmer kmer)

cdef class PSequenceFunction(PFunction):
    cpdef float evaluate(self, Sequence sequence)

