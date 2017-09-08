from numpy.random import ranf
from cython.operator cimport dereference as deref

cdef class ProbabilisticGraph:

    def __cinit__(self, Hashgraph graph, GraphFunction func):

        self.graph = graph

        self.K = graph.ksize()

        self.func = func
        self.func.set_graph(self.graph)


    cdef bool _insert_kmer(self, Kmer kmer):
        cdef float p = self.func.evaluate_kmer(kmer)
        if p > ranf():
            deref(self.graph._hg_this).add(deref(kmer._this).kmer_u)
            return True
        return False

    cdef bool _insert_sequence(self, Sequence sequence):
        cdef float p = self.func.evaluate_sequence(sequence)
        if p > ranf():
            deref(self.graph._hg_this).consume_string(sequence._obj.sequence)
            return True
        return False

    def insert(self, object item):
        if isinstance(item, Kmer):
            return self._insert_kmer(item)
        elif isinstance(item, Sequence):
            return self._insert_sequence(item)
        else:
            return False
