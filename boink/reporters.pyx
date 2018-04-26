# reporters.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref
from libcpp.memory cimport unique_ptr, make_unique

from boink.utils cimport _bstring

from boink.events cimport _EventListener


cdef class Reporter(EventListener):

    def __cinit__(self, str output_filename, *args, **kwargs):
        self.output_filename = output_filename
        if type(self) is Reporter:
            self._reporter = make_unique[_Reporter](_bstring(output_filename))
            self._listener = <_EventListener*>self._reporter.get()


cdef class StreamingCompactorReporter(Reporter):

    def __cinit__(self, str output_filename, *args, **kwargs):
        if type(self) is StreamingCompactorReporter:
            self._sc_reporter = make_unique[_StreamingCompactorReporter](_bstring(output_filename))
            self._reporter.reset(self._sc_reporter.get())
            self._listener = <_EventListener*>self._reporter.get()
