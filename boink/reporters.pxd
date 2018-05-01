# reporters.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libc.stdint cimport uint8_t, uint32_t, uint64_t
from libcpp.memory cimport shared_ptr, unique_ptr
from libcpp.string cimport string

from boink.dbg cimport DefaultDBG
from boink.cdbg cimport cDBGFormat, _cDBG, _StreamingCompactor
from boink.events cimport EventListener, _EventListener


cdef extern from "boink/reporters.hh" namespace "boink::reporters" nogil:
    cdef cppclass _SingleFileReporter "boink::reporters::SingleFileReporter" (_EventListener):
        _SingleFileReporter(string&)

    cdef cppclass _MultiFileReporter "boink::reporters::MultiFileReporter" (_EventListener):
        _MultiFileReporter(string&)

    cdef cppclass _StreamingCompactorReporter "boink::reporters::StreamingCompactorReporter" [GraphType] (_SingleFileReporter):
        _StreamingCompactorReporter(_StreamingCompactor[GraphType] *, string&)

    cdef cppclass _cDBGWriter "boink::reporters::cDBGWriter" (_MultiFileReporter):
        _cDBGWriter(_cDBG *, cDBGFormat, const string&)


cdef class SingleFileReporter(EventListener):
    cdef readonly object          output_filename
    cdef unique_ptr[_SingleFileReporter]    _owner
    cdef _SingleFileReporter *    _this


cdef class MultiFileReporter(EventListener):
    cdef readonly object                prefix
    cdef unique_ptr[_MultiFileReporter] _owner
    cdef _MultiFileReporter *           _this


cdef class StreamingCompactorReporter(SingleFileReporter):
    cdef unique_ptr[_StreamingCompactorReporter[DefaultDBG]] _s_owner
    cdef _StreamingCompactorReporter[DefaultDBG] *           _s_this


cdef class cDBGWriter(MultiFileReporter):
    cdef unique_ptr[_cDBGWriter] _s_owner
    cdef _cDBGWriter *           _s_this
