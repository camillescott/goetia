# reporters.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import sys

from cython.operator cimport dereference as deref
from libcpp.memory cimport unique_ptr, make_unique

from boink.cdbg cimport cDBG, StreamingCompactor, convert_format
from boink.utils cimport _bstring
from boink.events cimport _EventListener


cdef class Reporter(EventListener):

    def __cinit__(self, str output_filename,
                        *args, **kwargs):
        self.output_filename = output_filename
        if type(self) is Reporter:
            self._reporter_owner = make_unique[_Reporter](_bstring(output_filename))
            self._reporter_this = self._reporter_owner.get()
            self._listener = <_EventListener*>self._reporter_owner.get()


cdef class StreamingCompactorReporter(Reporter):

    def __cinit__(self, str output_filename, StreamingCompactor compactor,
                        *args, **kwargs):
        if type(self) is StreamingCompactorReporter:
            self._sc_reporter_owner = make_unique[_StreamingCompactorReporter[DefaultDBG]](\
                    compactor._sc_this.get(), _bstring(output_filename))
            self._sc_reporter_this = self._sc_reporter_owner.get()
            self._reporter_this    = self._sc_reporter_this
            self._listener = <_EventListener*>self._sc_reporter_owner.get()

    def __dealloc__(self):
        print('Dealloc StreamingCompactorReporter', file=sys.stderr)


cdef class cDBGWriter(Reporter):

    def __cinit__(self, str output_filename, str graph_format, cDBG cdbg,
                        *args, **kwargs):
        if type(self) is cDBGWriter:
            self._writer_reporter_owner = make_unique[_cDBGWriter](cdbg._this,
                                                             convert_format(graph_format),
                                                             _bstring(output_filename))
            self._writer_reporter_this = self._writer_reporter_owner.get()
            self._reporter_this = self._writer_reporter_this
            self._listener = <_EventListener*>self._writer_reporter_owner.get()
        
