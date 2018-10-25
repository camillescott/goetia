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
from boink.cdbg cimport cDBGFormat, _cDBG 
from boink.compactor cimport *
from boink.events cimport EventListener, _EventListener


cdef extern from "boink/reporters.hh" namespace "boink::reporters" nogil:
    cdef cppclass _SingleFileReporter "boink::reporters::SingleFileReporter" (_EventListener):
        _SingleFileReporter(string&)

    cdef cppclass _MultiFileReporter "boink::reporters::MultiFileReporter" (_EventListener):
        _MultiFileReporter(string&)

    cdef cppclass _StreamingCompactorReporter "boink::reporters::StreamingCompactorReporter" [GraphType] (_SingleFileReporter):
        _StreamingCompactorReporter(_StreamingCompactor[GraphType] *, string&)

    cdef cppclass _cDBGWriter "boink::reporters::cDBGWriter" [GraphType] (_MultiFileReporter):
        _cDBGWriter(_cDBG[GraphType] *, cDBGFormat, const string&)

    cdef cppclass _cDBGHistoryReporter "boink::reporters::cDBGHistoryReporter" (_SingleFileReporter):
        _cDBGHistoryReporter(const string&)


cdef class SingleFileReporter(EventListener):
    cdef readonly object          output_filename
    cdef unique_ptr[_SingleFileReporter]    _owner
    cdef _SingleFileReporter *    _this


cdef class MultiFileReporter(EventListener):
    cdef readonly object                prefix
    cdef unique_ptr[_MultiFileReporter] _owner
    cdef _MultiFileReporter *           _this


cdef class cDBGHistoryReporter(SingleFileReporter):
    cdef unique_ptr[_cDBGHistoryReporter] _h_owner
    cdef _cDBGHistoryReporter *           _h_this

include "reporters.tpl.pxd.pxi"
