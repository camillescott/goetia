# cython: c_string_type=unicode, c_string_encoding=utf8
import cython
from cython.operator cimport dereference as deref, preincrement as inc

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.set cimport set
from libcpp.queue cimport queue
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr, make_shared

from libc.stdint cimport uint8_t, uint32_t, uint64_t

from khmer._oxli.partitioning cimport StreamingPartitioner, Component
from khmer._oxli.parsing cimport Sequence
from khmer._oxli.wrapper cimport *

from boink.stats cimport PartitionFunction


cdef class ConditionalPartitioner(StreamingPartitioner):

    def __cinit__(self, graph, tag_density=None):
        self.min_fitness = 0.95

    def consume(self, Sequence sequence, PartitionFunction func=None):
        if func is None:
            return deref(self._this).consume(sequence._obj.sequence)
        cdef float fitness
        return self._consume_conditional(sequence._obj.sequence, func, fitness)
        
    cdef int _consume_conditional(self, string sequence, PartitionFunction func, 
                                  float& fitness) except -1:

        cdef set[HashIntoType] tags
        cdef KmerQueue seeds
        cdef set[HashIntoType] seen
        cdef uint64_t n_new = deref(self._this).seed_sequence(sequence, tags, seeds, seen)
        (&fitness)[0] = func._evaluate_tags(sequence, tags)
        if fitness > self.min_fitness: # blerp
            deref(self._this).find_connected_tags(seeds, tags, seen, False)
            deref(self._this).create_and_connect_components(tags)

        return n_new
