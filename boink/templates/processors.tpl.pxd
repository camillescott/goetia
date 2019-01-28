{# boink/templates/processors.tpl.pxd
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}
{% extends "base.tpl" %}
{% block code %}

from libc.stdint cimport uint16_t, uint32_t, uint64_t, int64_t
from libcpp cimport bool
from libcpp.string cimport string

from boink.dbg cimport *
from boink.cdbg cimport *
from boink.events cimport EventNotifier, _EventNotifier, _EventListener
from libcpp.memory cimport shared_ptr
from boink.utils cimport _bstring

from boink.processors cimport *

cdef class FileProcessor_Base(FileProcessor):
    cdef readonly object storage_type
    cdef readonly object shifter_type

{% for type_bundle in type_bundles %}

cdef class FileConsumer_{{type_bundle.suffix}}(FileProcessor_Base):
    cdef shared_ptr[_FileConsumer[_dBG[{{type_bundle.params}}]]] _this


cdef class DecisionNodeProcessor_{{type_bundle.suffix}}(FileProcessor_Base):
    cdef readonly str output_filename
    cdef shared_ptr[_DecisionNodeProcessor[_dBG[{{type_bundle.params}}]]] _this


cdef class StreamingCompactorProcessor_{{type_bundle.suffix}}(FileProcessor_Base):
    cdef readonly str output_filename
    cdef shared_ptr[_StreamingCompactorProcessor[_dBG[{{type_bundle.params}}]]] _this


cdef class NormalizingCompactor_{{type_bundle.suffix}}(FileProcessor_Base):
    cdef readonly str output_filename
    cdef shared_ptr[_NormalizingCompactor[_dBG[{{type_bundle.params}}]]] _this


{% endfor %}

cdef object _make_file_consumer(dBG_Base graph,
                                uint64_t fine_interval,
                                uint64_t medium_interval,
                                uint64_t coarse_interval)

cdef object _make_decision_node_processor(StreamingCompactor_Base compactor,
                                          str filename, 
                                          uint64_t fine_interval,
                                          uint64_t medium_interval,
                                          uint64_t coarse_interval)

cdef object _make_streaming_compactor_processor(StreamingCompactor_Base compactor, 
                                                uint64_t fine_interval,
                                                uint64_t medium_interval,
                                                uint64_t coarse_interval)

cdef object _make_normalizing_compactor(StreamingCompactor_Base compactor, 
                                        unsigned int cutoff,
                                        uint64_t fine_interval,
                                        uint64_t medium_interval,
                                        uint64_t coarse_interval)


{% endblock code %}
