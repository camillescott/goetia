from libcpp.memory cimport shared_ptr, weak_ptr

from khmer._oxli.wrapper cimport CpHashgraph
from khmer._oxli.hashing cimport Kmer
from khmer._oxli.parsing cimport Sequence
from khmer._oxli.traversal cimport Traverser

from cython cimport numeric


cdef class GraphFunction:
    # have PFunction store Hashgraph ptr
    cdef CpHashgraph * graph
    cdef void _set_graph(self, CpHashgraph *)
    cpdef float evaluate_kmer(self, Kmer)
    cpdef float evaluate_sequence(self, Sequence)

cdef class KmerCountFunction(GraphFunction):
    cpdef float evaluate_kmer(self, Kmer)

cdef class KmerDegreeFunction(GraphFunction):
    cdef Traverser traverser
    cpdef float evaluate_kmer(self, Kmer)

