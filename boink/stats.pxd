cimport cython

from libcpp.memory cimport shared_ptr, weak_ptr
from libcpp.string cimport string
from libcpp.set cimport set
from libcpp.vector cimport vector

from khmer._oxli.oxli_types cimport *
from khmer._oxli.graphs cimport CpHashgraph, Hashgraph
from khmer._oxli.partitioning cimport CpStreamingPartitioner
from khmer._oxli.assembly cimport CompactingAssembler
from khmer._oxli.hashing cimport Kmer, CpKmer
from khmer._oxli.sequence cimport Sequence
from khmer._oxli.traversal cimport Traverser
from libc.stdint cimport uint8_t, uint16_t, uint64_t


cdef extern from 'constants.h':
    enum: N_HEPTAMERS
    enum: MAX_DBG_DEGREE


cdef class GraphFunction:
    # have PFunction store Hashgraph ptr
    cdef Hashgraph graph
    cdef public WordLength K
    cpdef int consume(self, Sequence sequence)
    cpdef float evaluate_kmer(self, Kmer)
    cpdef float evaluate_sequence(self, Sequence)


cdef class PartitionFunction(GraphFunction):
    cdef shared_ptr[CpStreamingPartitioner] partitioner
    cdef float _evaluate_tags(self, string,
                              vector[HashIntoType]&)
    cpdef float evaluate_tags(self, Sequence, tuple)


cdef class PartitionCoverage(PartitionFunction):
    cdef float coverage_cutoff


cdef class PartitionCoverageSlice(PartitionFunction):
    cdef float coverage_floor
    cdef float coverage_ceiling


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
