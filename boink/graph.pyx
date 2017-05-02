from khmer._oxli.wrapper cimport get_hashgraph_ptr
from khmer._oxli.hashing import Kmer
from khmer._oxli.parsing import Sequence

from numpy.random import ranf
from cython.operator cimport dereference as deref

cdef class ProbabilisticGraph:

    def __cinit__(self, object graph, GraphFunction func):

        self._graph = get_hashgraph_ptr(graph)

        self.K = deref(self._graph).ksize()

        self.func = func
        self.func._set_graph(self._graph)


    cdef bool _insert_kmer(self, Kmer kmer):
        cdef float p = self.func.evaluate_kmer(kmer)
        if p > ranf():
            deref(self._graph).add(deref(kmer._this).kmer_u)
            return True
        return False

    cdef bool _insert_sequence(self, Sequence sequence):
        cdef float p = self.func.evaluate_sequence(sequence)
        if p > ranf():
            deref(self._graph).consume_string(sequence._obj.sequence)
            return True
        return False

    def insert(self, object item):
        if isinstance(item, Kmer):
            return self._insert_kmer(item)
        elif isinstance(item, Sequence):
            return self._insert_sequence(item)
        else:
            return False
