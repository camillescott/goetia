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
    cdef cppclass _Reporter "boink::reporters::Reporter" (_EventListener):
        _Reporter(string&)

    cdef cppclass _StreamingCompactorReporter "boink::reporters::StreamingCompactorReporter" [GraphType] (_Reporter):
        _StreamingCompactorReporter(_StreamingCompactor[GraphType] *, string&)

    cdef cppclass _cDBGWriter "boink::reporters::cDBGWriter" (_Reporter):
        _cDBGWriter(_cDBG *, cDBGFormat, const string&)


cdef class Reporter(EventListener):
    cdef readonly object       output_filename
    cdef unique_ptr[_Reporter] _owner
    cdef _Reporter *           _this


cdef class StreamingCompactorReporter(Reporter):
    cdef unique_ptr[_StreamingCompactorReporter[DefaultDBG]] _s_owner
    cdef _StreamingCompactorReporter[DefaultDBG] *           _s_this


cdef class cDBGWriter(Reporter):
    cdef unique_ptr[_cDBGWriter] _s_owner
    cdef _cDBGWriter *           _s_this
