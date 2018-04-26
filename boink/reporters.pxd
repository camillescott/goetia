# reporters.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libc.stdint cimport uint8_t, uint32_t, uint64_t
from libcpp.memory cimport shared_ptr, unique_ptr
from libcpp.string cimport string

from boink.events cimport EventListener, _EventListener


cdef extern from "boink/reporters.hh" namespace "boink::reporters" nogil:
    cdef cppclass _Reporter "boink::reporters::Reporter" (_EventListener):
        _Reporter(string&)

    cdef cppclass _StreamingCompactorReporter "boink::reporters::StreamingCompactorReporter" (_Reporter):
        _StreamingCompactorReporter(string&)


cdef class Reporter(EventListener):
    cdef readonly object output_filename
    cdef unique_ptr[_Reporter] _reporter


cdef class StreamingCompactorReporter(Reporter):
    cdef unique_ptr[_StreamingCompactorReporter] _sc_reporter
