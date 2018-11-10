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


class DEFAULT_INTERVALS:
    FINE   = DEFAULT_FINE_INTERVAL
    MEDIUM = DEFAULT_MEDIUM_INTERVAL
    COARSE = DEFAULT_COARSE_INTERVAL


cdef class FileProcessor:

    def __cinit__(self, *args, **kwargs):
        pass

include "processors.tpl.pyx.pxi"

cdef class MinimizerProcessor(FileProcessor):

    def __cinit__(self, int64_t window_size,
                        uint16_t ksize,
                        str output_filename,
                        uint64_t fine_interval,
                        uint64_t medium_interval,
                        uint64_t coarse_interval):

        cdef string _output_filename = _bstring(output_filename)
        self._mp_this = make_unique[_MinimizerProcessor[_DefaultShifter]](window_size,
                                                                         ksize,
                                                                         _output_filename,
                                                                         fine_interval,
                                                                         medium_interval,
                                                                         coarse_interval)
    def process(self, str input_filename):
        deref(self._mp_this).process(_bstring(input_filename))

        return deref(self._mp_this).n_reads()


def make_file_consumer(dBG_Base graph,
                       int fine_interval=DEFAULT_FINE_INTERVAL,
                       int medium_interval=DEFAULT_MEDIUM_INTERVAL,
                       int coarse_interval=DEFAULT_COARSE_INTERVAL):

    return _make_file_consumer(graph,
                               fine_interval,
                               medium_interval,
                               coarse_interval)


def make_decision_node_processor(StreamingCompactor_Base compactor,
                                 str output_filename,
                                 int fine_interval=DEFAULT_FINE_INTERVAL,
                                 int medium_interval=DEFAULT_MEDIUM_INTERVAL,
                                 int coarse_interval=DEFAULT_COARSE_INTERVAL):

    return _make_decision_node_processor(compactor,
                                         output_filename,
                                         fine_interval,
                                         medium_interval,
                                         coarse_interval)


def make_streaming_compactor_processor(StreamingCompactor_Base compactor,
                                       int fine_interval=DEFAULT_FINE_INTERVAL,
                                       int medium_interval=DEFAULT_MEDIUM_INTERVAL,
                                       int coarse_interval=DEFAULT_COARSE_INTERVAL):
        
    return _make_streaming_compactor_processor(compactor,
                                               fine_interval,
                                               medium_interval,
                                               coarse_interval)
