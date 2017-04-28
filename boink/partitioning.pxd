# cython: c_string_type=unicode, c_string_encoding=utf8
import cython
from cython.operator cimport dereference as deref, preincrement as inc
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.set cimport set
from libcpp.map cimport map
from libcpp.queue cimport queue
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr, make_shared

from libc.stdint cimport uint8_t, uint32_t, uint64_t

from khmer._oxli.wrapper cimport HashIntoType, KmerQueue, CpKmer, CpHashgraph
from khmer._oxli.wrapper cimport CpStreamingPartitioner
from khmer._oxli.partitioning cimport StreamingPartitioner
from boink.stats cimport PartitionFunction

__CythonIsBroken__ = MemoryError


cdef class ConditionalPartitioner(StreamingPartitioner):

    cdef float min_fitness
    cdef int _consume_conditional(self, string, PartitionFunction, float&) except -1
    cdef int _consume_conditional_pair(self, string, string,
                                       PartitionFunction, float&) except -1

cdef class DoConditionalPartitioning:
    cdef object args
    cdef readonly object graph
    cdef readonly ConditionalPartitioner partitioner


