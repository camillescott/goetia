from libcpp.memory cimport shared_ptr, weak_ptr
from cython cimport numeric
from cython.operator cimport dereference as deref

from khmer._oxli.wrapper cimport CpHashgraph, get_hashgraph_ptr
from khmer._oxli.hashing cimport Kmer
from khmer._oxli.hashing import Kmer
from khmer._oxli.parsing cimport Sequence
from khmer._oxli.parsing import Sequence
from khmer._oxli.traversal cimport Traverser


cdef class GraphFunction:

    def __cinit__(self, object graph=None, *args, **kwargs):
        if graph is not None:
            self.graph = get_hashgraph_ptr(graph)

    cdef void _set_graph(self, CpHashgraph * ptr):
        self.graph = ptr

    def set_graph(self, object graph):
        self.graph = get_hashgraph_ptr(graph)

    cpdef float evaluate_kmer(self, Kmer kmer):
        raise NotImplementedError()

    cpdef float evaluate_sequence(self, Sequence sequence):
        raise NotImplementedError()

    def save(self, filename):
        raise NotImplementedError()

    @staticmethod
    def load(filename):
        raise NotImplementedError()


cdef class KmerCountFunction(GraphFunction):

    cpdef float evaluate_kmer(self, Kmer kmer):
        return <int>deref(self.graph).get_count(deref(kmer._this.get()).kmer_u)


cdef class KmerDegreeFunction(GraphFunction):

    def __cinit__(self, object graph=None, *args, **kwargs):
        self.traverser = Traverser(graph)

    def set_graph(self, object graph):
        self.graph = get_hashgraph_ptr(graph)
        self.traverser = Traverser(graph)

    cpdef float evaluate_kmer(self, Kmer kmer):
        return <int>self.traverser.degree(kmer)


def GraphFunction_shim(GraphFunction func, object graph, *args, **kwargs):
    cdef CpHashgraph * ptr = get_hashgraph_ptr(graph)
    func._set_graph(ptr)
    return func.evaluate(*args, **kwargs)
