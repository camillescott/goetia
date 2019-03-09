# reporters.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import sys

from cython.operator cimport dereference as deref
from libcpp.memory cimport shared_ptr, make_shared

from boink.cdbg cimport *
from boink.utils cimport _bstring
from boink.events cimport _EventListener
from boink.minimizers cimport UKHSCountSignature


cdef class SingleFileReporter(EventListener):

    def __cinit__(self, str output_filename,
                        *args, **kwargs):
        self.output_filename = output_filename
        if type(self) is SingleFileReporter:
            self._this = make_shared[_SingleFileReporter](_bstring(output_filename))
            self._listener = <shared_ptr[_EventListener]>self._this


cdef class MultiFileReporter(EventListener):

    def __cinit__(self, str prefix,
                        *args, **kwargs):
        self.prefix = prefix
        if type(self) is MultiFileReporter:
            self._this = make_shared[_MultiFileReporter](_bstring(prefix))
            self._listener = <shared_ptr[_EventListener]>self._this


cdef class cDBGHistoryReporter(SingleFileReporter):

    def __cinit__(self, str output_filename, *args, **kwargs):
        if type(self) is cDBGHistoryReporter:
            self._h_this = make_shared[_cDBGHistoryReporter](_bstring(output_filename));
            self._this = <shared_ptr[_SingleFileReporter]>self._h_this
            self._listener = <shared_ptr[_EventListener]>self._h_this

cdef class UKHSSignatureReporter(SingleFileReporter):


    def __cinit__(self, str output_filename, UKHSCountSignature signature):
        if type(self) is UKHSSignatureReporter:
            self._uk_this = make_shared[_UKHSSignatureReporter](signature._this,
                                                                _bstring(output_filename))
            self._this = <shared_ptr[_SingleFileReporter]>self._uk_this
            self._listener = <shared_ptr[_EventListener]>self._uk_this

include "reporting.tpl.pyx.pxi"

