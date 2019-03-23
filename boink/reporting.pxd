# reporters.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libc.stdint cimport uint8_t, uint32_t, uint64_t
from libcpp.memory cimport shared_ptr
from libcpp.string cimport string

from boink.cdbg cimport cDBGFormat, _cDBG 
from boink.compactor cimport *
from boink.events cimport EventListener, _EventListener
from boink.minimizers cimport _UKHSCountSignature
from boink.prometheus cimport _Registry


cdef extern from "boink/reporting/reporters.hh" namespace "boink::reporting" nogil:
    cdef cppclass _SingleFileReporter "boink::reporting::SingleFileReporter" (_EventListener):
        _SingleFileReporter(string&)

    cdef cppclass _MultiFileReporter "boink::reporting::MultiFileReporter" (_EventListener):
        _MultiFileReporter(string&)

cdef extern from "boink/reporting/streaming_compactor_reporter.hh" namespace "boink::reporting" nogil:
    cdef cppclass _StreamingCompactorReporter "boink::reporting::StreamingCompactorReporter" [GraphType] (_SingleFileReporter):
        _StreamingCompactorReporter(shared_ptr[_StreamingCompactor[GraphType]], string&)

cdef extern from "boink/reporting/cdbg_writer_reporter.hh" namespace "boink::reporting" nogil:
    cdef cppclass _cDBGWriter "boink::reporting::cDBGWriter" [GraphType] (_MultiFileReporter):
        _cDBGWriter(shared_ptr[_cDBG[GraphType]], cDBGFormat, const string&)

cdef extern from "boink/reporting/cdbg_history_reporter.hh" namespace "boink::reporting" nogil:
    cdef cppclass _cDBGHistoryReporter "boink::reporting::cDBGHistoryReporter" (_SingleFileReporter):
        _cDBGHistoryReporter(const string&)

cdef extern from "boink/reporting/ukhs_signature_reporter.hh" namespace "boink::reporting" nogil:
    cdef cppclass _UKHSSignatureReporter "boink::reporting::UKHSSignatureReporter" (_SingleFileReporter):
        _UKHSSignatureReporter(shared_ptr[_UKHSCountSignature], const string&)

cdef extern from "boink/reporting/cdbg_component_reporter.hh" namespace "boink::reporting" nogil:
    cdef cppclass _cDBGComponentReporter "boink::reporting::cDBGComponentReporter" [GraphType] (_SingleFileReporter):
        _cDBGComponentReporter(shared_ptr[_cDBG[GraphType]],
                               const string&,
                               shared_ptr[_Registry])
        _cDBGComponentReporter(shared_ptr[_cDBG[GraphType]],
                               const string&,
                               shared_ptr[_Registry]
                               size_t)

cdef extern from "boink/reporting/cdbg_unitig_reporter.hh" namespace "boink::reporting" nogil:
    cdef cppclass _cDBGUnitigReporter "boink::reporting::cDBGUnitigReporter" [GraphType] (_SingleFileReporter):
        _cDBGUnitigReporter(shared_ptr[_cDBG[GraphType]],
                            const string&,
                            vector[size_t])


cdef class SingleFileReporter(EventListener):
    cdef readonly object          output_filename
    cdef shared_ptr[_SingleFileReporter]    _this


cdef class MultiFileReporter(EventListener):
    cdef readonly object                prefix
    cdef shared_ptr[_MultiFileReporter] _this


cdef class cDBGHistoryReporter(SingleFileReporter):
    cdef shared_ptr[_cDBGHistoryReporter] _h_this


cdef class UKHSSignatureReporter(SingleFileReporter):
    cdef shared_ptr[_UKHSSignatureReporter] _uk_this

include "reporting.tpl.pxd.pxi"
