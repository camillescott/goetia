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


cdef class FileConsumer(FileProcessor):

    def __cinit__(self, dBG__BitStorage__DefaultShifter graph,
                        uint32_t output_interval):
        self._fc_this = make_unique[_FileConsumer[DefaultDBG]](graph._this.get(),
                                                               output_interval)

    def process(self, str input_filename):
        deref(self._fc_this).process(_bstring(input_filename))

        return (deref(self._fc_this).n_reads(),
                deref(self._fc_this).n_consumed())



cdef class DecisionNodeProcessor(FileProcessor):

    def __cinit__(self, StreamingCompactor compactor, str output_filename,
                        uint32_t output_interval):
        self.output_filename = output_filename
        cdef string _output_filename = _bstring(output_filename)
        self._dnp_this = make_unique[_DecisionNodeProcessor[DefaultDBG]](compactor._sc_this.get(),
                                                                         _output_filename,
                                                                         output_interval)

    def process(self, str input_filename):
        deref(self._dnp_this).process(_bstring(input_filename))

        return deref(self._dnp_this).n_reads()


cdef class StreamingCompactorProcessor(FileProcessor):

    def __cinit__(self, StreamingCompactor compactor, str output_filename,
                        uint32_t output_interval):
        self.output_filename = output_filename
        cdef string _output_filename = _bstring(output_filename)
        self._scp_this = make_unique[_StreamingCompactorProcessor[DefaultDBG]](compactor._sc_this.get(),
                                                                               _output_filename,
                                                                               output_interval)

    def process(self, str input_filename):
        deref(self._scp_this).process(_bstring(input_filename))

        return deref(self._scp_this).n_reads()


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



