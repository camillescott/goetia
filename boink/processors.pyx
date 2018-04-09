# processors.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref
from libcpp.memory cimport make_unique

from boink.dbg cimport *
from boink.cdbg cimport *
from boink.utils cimport _bstring, _ustring


cdef class FileProcessor:

    def __cinit__(self, *args, **kwargs):
        pass


cdef class FileConsumer:

    def __cinit__(self, dBG_BitStorage_DefaultShifter graph):
        self._fc_this = make_unique[_FileConsumer[DefaultDBG]](graph._this.get())

    def process(self, str input_filename):
        cdef uint64_t n_reads = 0
        cdef uint64_t n_kmers = 0

        deref(self._fc_this).process(_bstring(input_filename),
                                     n_reads,
                                     n_kmers)

        return n_reads, n_kmers


cdef class DecisionNodeProcessor:

    def __cinit__(self, StreamingCompactor compactor, str output_filename):
        self.output_filename = output_filename
        cdef string _output_filename = _bstring(output_filename)
        self._dnp_this = make_unique[_DecisionNodeProcessor[DefaultDBG]](compactor._sc_this.get(),
                                                                         _output_filename)

    def process(self, str input_filename):
        cdef uint64_t n_reads = 0
        cdef uint64_t n_kmers = 0

        deref(self._dnp_this).process(_bstring(input_filename),
                                      n_reads,
                                      n_kmers)

        return n_reads, n_kmers
