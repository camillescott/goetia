{# boink/templates/reporters.tpl.pxd
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}
{% extends "base.tpl" %}
{% block code %}

from libc.stdint cimport uint8_t, uint32_t, uint64_t
from libcpp.memory cimport shared_ptr
from libcpp.string cimport string

from boink.dbg cimport *
from boink.compactor cimport *
from boink.cdbg cimport cDBGFormat, _cDBG
from boink.events cimport EventListener, _EventListener

from boink.processors cimport *


cdef class StreamingCompactorReporter_Base(SingleFileReporter):
    cdef readonly object storage_type
    cdef readonly object shifter_type


cdef class cDBGWriter_Base(MultiFileReporter):
    cdef readonly object storage_type
    cdef readonly object shifter_type


cdef class cDBGComponentReporter(SingleFileReporter):
    cdef readonly object storage_type
    cdef readonly object shifter_type


cdef class cDBGUnitigReporter(SingleFileReporter):
    cdef readonly object storage_type
    cdef readonly object shifter_type


{% for type_bundle in type_bundles %}
cdef class StreamingCompactorReporter_{{type_bundle.suffix}}(StreamingCompactorReporter_Base):
    cdef shared_ptr[_StreamingCompactorReporter[_dBG[{{type_bundle.params}}]]] _s_this

cdef class cDBGWriter_{{type_bundle.suffix}}(cDBGWriter_Base):
    cdef shared_ptr[_cDBGWriter[_dBG[{{type_bundle.params}}]]] _s_this

cdef class cDBGComponentReporter_{{type_bundle.suffix}}(cDBGComponentReporter):
    cdef shared_ptr[_cDBGComponentReporter[_dBG[{{type_bundle.params}}]]] _s_this

cdef class cDBGUnitigReporter_{{type_bundle.suffix}}(cDBGUnitigReporter):
    cdef shared_ptr[_cDBGUnitigReporter[_dBG[{{type_bundle.params}}]]] _s_this


{% endfor %}

cdef object _make_streaming_compactor_reporter(str output_filename,
                                               StreamingCompactor compactor)

cdef object _make_cdbgwriter_reporter(str output_prefix,
                                      str graph_format,
                                      cDBG_Base cdbg)

{% endblock code %}
