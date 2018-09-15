# reporters.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import sys

from cython.operator cimport dereference as deref
from libcpp.memory cimport unique_ptr, make_unique

from boink.cdbg cimport cDBG, convert_format
from boink.utils cimport _bstring
from boink.events cimport _EventListener


cdef class SingleFileReporter(EventListener):

    def __cinit__(self, str output_filename,
                        *args, **kwargs):
        self.output_filename = output_filename
        if type(self) is SingleFileReporter:
            self._owner = make_unique[_SingleFileReporter](_bstring(output_filename))
            self._this = self._owner.get()
            self._listener = <_EventListener*>self._owner.get()


cdef class MultiFileReporter(EventListener):

    def __cinit__(self, str prefix,
                        *args, **kwargs):
        self.prefix = prefix
        if type(self) is MultiFileReporter:
            self._owner = make_unique[_MultiFileReporter](_bstring(prefix))
            self._this = self._owner.get()
            self._listener = <_EventListener*>self._owner.get()


cdef class cDBGWriter(MultiFileReporter):

    def __cinit__(self, str output_prefix, str graph_format, cDBG cdbg,
                        *args, **kwargs):
        if type(self) is cDBGWriter:
            self._s_owner = make_unique[_cDBGWriter](cdbg._this,
                                                     convert_format(graph_format),
                                                     _bstring(output_prefix))
            self._s_this = self._s_owner.get()
            self._this = self._s_this
            self._listener = <_EventListener*>self._s_owner.get()

include "reporters.tpl.pyx.pxi"

def make_streaming_compactor_reporter(str output_filename, StreamingCompactor_Base compactor):
    return _make_streaming_compactor_reporter(output_filename, compactor)
