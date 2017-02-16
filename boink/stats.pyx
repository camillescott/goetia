from libcpp.memory cimport shared_ptr, weak_ptr

from khmer._oxli.wrapper cimport CpHashgraph, get_hashgraph_ptr
from khmer._oxli.hashing cimport Kmer
from khmer._oxli.hashing import Kmer
from khmer._oxli.parsing cimport Sequence
from khmer._oxli.parsing import Sequence


cdef class PFunction:

    def __cinit__(self, object graph=None):
        if graph is not None:
            self.graph = get_hashgraph_ptr(graph)

    cdef void _set_graph(self, CpHashgraph * ptr):
        self.graph = ptr

    def set_graph(self, object graph):
        self.graph = get_hashgraph_ptr(graph)

    def save(self, filename):
        raise NotImplementedError()

    @staticmethod
    def load(filename):
        raise NotImplementedError()


cdef class PKmerFunction(PFunction):

    cpdef float evaluate(self, Kmer kmer):
        return 1.0


cdef class PSequenceFunction(PFunction):

    cpdef float evaluate(self, Sequence sequence):
        return 1.0






def PFunction_shim(PFunction func, object graph, *args, **kwargs):
    cdef CpHashgraph * ptr = get_hashgraph_ptr(graph)
    func._set_graph(ptr)
    return func.evaluate(*args, **kwargs)
