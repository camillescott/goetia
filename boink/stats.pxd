from libcpp.memory cimport shared_ptr, weak_ptr
from libcpp.string cimport string
cimport cython

from khmer._oxli.wrapper cimport CpHashgraph, HashIntoType, WordLength, CpStreamingPartitioner
from khmer._oxli.assembly cimport CompactingAssembler
from khmer._oxli.hashing cimport Kmer
from khmer._oxli.parsing cimport Sequence
from khmer._oxli.traversal cimport Traverser
from libc.stdint cimport uint8_t, uint16_t, uint64_t
from libcpp.set cimport set


cdef extern from 'constants.h':
    enum: N_HEPTAMERS
    enum: MAX_DBG_DEGREE


cdef class GraphFunction:
    # have PFunction store Hashgraph ptr
    cdef CpHashgraph * graph
    cdef public WordLength K
    cpdef int consume(self, Sequence sequence)
    cdef void _set_graph(self, CpHashgraph *)
    cpdef float evaluate_kmer(self, Kmer)
    cpdef float evaluate_sequence(self, Sequence)


cdef class PartitionFunction(GraphFunction):
    cdef shared_ptr[CpStreamingPartitioner] partitioner
    cdef float _evaluate_tags(self, string,
                              set[HashIntoType]&)
    cpdef float evaluate_tags(self, Sequence, tuple)


cdef class PartitionCoverage(PartitionFunction):
    cdef float coverage_cutoff


cdef class KmerCountFunction(GraphFunction):
    cpdef float evaluate_kmer(self, Kmer)


cdef class KmerDegreeFunction(GraphFunction):
    cdef Traverser traverser


cdef class KmerDegreeBiasFunction(KmerDegreeFunction):

    cdef uint64_t heptamer_counts[N_HEPTAMERS][MAX_DBG_DEGREE]
    cdef uint64_t degree_counts[MAX_DBG_DEGREE]
    #cpdef float evaluate_kmer(self, Kmer)


cdef class SamplePathsFunction(GraphFunction):

    cdef CompactingAssembler assembler
    cdef uint64_t bandwidth
    cdef uint64_t n_batches
    cdef list results
    cdef uint64_t t
