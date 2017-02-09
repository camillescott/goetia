from khmer._oxli.utils cimport _get_n_primes_near_x
from khmer._oxli._oxli cimport get_hashgraph_ptr
from khmer._oxli.hashing import Kmer
from khmer._oxli.parsing import Sequence

from numpy.random import ranf
from cython.operator cimport dereference as deref

cdef class ProbabilisticGraph:

    def __cinit__(self, object graph,
                        PKmerFunction kmer_func, 
                        PSequenceFunction sequence_func):

        self._graph = get_hashgraph_ptr(graph)

        self.K = deref(self._graph).ksize()
        self.kmer_func = kmer_func
        self.sequence_func = sequence_func


    cdef bool _insert_kmer(self, Kmer kmer):
        cdef float p = self.kmer_func.evaluate(kmer, self._graph)
        if p < ranf():
            deref(self._graph).add(deref(kmer._this).kmer_u)
            return True
        return False

    cdef bool _insert_sequence(self, Sequence sequence):
        cdef float p = self.sequence_func.evaluate(sequence, self._graph)
        if p < ranf():
            deref(self._graph).consume_string(sequence._obj.sequence)
            return True
        return False


