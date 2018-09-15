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

include "processors.tpl.pyx.pxi"

cdef class MinimizerProcessor(FileProcessor):

    def __cinit__(self, int64_t window_size,
                        uint16_t ksize,
                        str output_filename,
                        uint32_t output_interval):
        cdef string _output_filename = _bstring(output_filename)
        self._mp_this = make_unique[_MinimizerProcessor[_DefaultShifter]](window_size,
                                                                         ksize,
                                                                         _output_filename,
                                                                         output_interval)
    def process(self, str input_filename):
        deref(self._mp_this).process(_bstring(input_filename))

        return deref(self._mp_this).n_reads()


def make_file_consumer(dBG_Base graph, int output_interval):
    return _make_file_consumer(graph, output_interval)


def make_decision_node_processor(StreamingCompactor_Base compactor, str output_filename, int output_interval):
    return _make_decision_node_processor(compactor, output_filename, output_interval)


def make_streaming_compactor_processor(StreamingCompactor_Base compactor, int output_interval):
    return _make_streaming_compactor_processor(compactor, output_interval)
