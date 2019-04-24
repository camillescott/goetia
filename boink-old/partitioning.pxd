# boink/partitioning.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref, preincrement as inc
from libcpp cimport bool
from libcpp.string cimport string
from libc.stdint cimport uint8_t, uint32_t, uint64_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.hashing cimport *
from khmer._oxli.graphs cimport CpHashgraph
from khmer._oxli.partitioning cimport CpStreamingPartitioner, StreamingPartitioner

from boink.stats cimport PartitionFunction

__CythonIsBroken__ = MemoryError


cdef class ConditionalPartitioner(StreamingPartitioner):

    cdef float min_fitness
    cdef readonly object info_filename
    cdef readonly object info_fp
    cdef int _consume_conditional(self, string, PartitionFunction, float&) except -1
    cdef int _consume_conditional_pair(self, string, string,
                                       PartitionFunction, float&) except -1

cdef class DoConditionalPartitioning:
    cdef object args
    cdef readonly object graph
    cdef readonly ConditionalPartitioner partitioner


