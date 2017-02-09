from libcpp.memory cimport shared_ptr, weak_ptr

from khmer._oxli._oxli cimport CpHashgraph
from khmer._oxli.hashing cimport Kmer
from khmer._oxli.hashing import Kmer
from khmer._oxli.parsing cimport Sequence
from khmer._oxli.parsing import Sequence


cdef class PKmerFunction:

    cdef float evaluate(self, Kmer kmer, CpHashgraph * graph):
        return 1.0


cdef class PSequenceFunction:

    cdef float evaluate(self, Sequence sequence, CpHashgraph * graph):
        return 1.0


