{# boink/templates/compactor.tpl.pxd
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}
{% extends "base.tpl" %}
{% block code %}

from libcpp.memory cimport shared_ptr

from boink.dbg cimport *
from boink.cdbg cimport *
from boink.prometheus cimport Instrumentation


cdef class StreamingCompactor:
    cdef readonly object storage_type
    cdef readonly object shifter_type
    cdef object graph

cdef class SolidStreamingCompactor:
    cdef readonly object    storage_type
    cdef readonly object    shifter_type
    cdef StreamingCompactor compactor
    cdef dBG                abund_filter

{% for type_bundle in type_bundles %}
cdef class StreamingCompactor_{{type_bundle.suffix}}(StreamingCompactor):
    cdef shared_ptr[_StreamingCompactor[_dBG[{{type_bundle.params}}]]] _this
    cdef shared_ptr[_dBG[{{type_bundle.params}}]] _graph
    cdef public cDBG_{{type_bundle.suffix}} cdbg
    cdef public EventNotifier Notifier
    cdef Instrumentation instrumentation

cdef class SolidStreamingCompactor_{{type_bundle.suffix}}(SolidStreamingCompactor):
    cdef shared_ptr[_SolidStreamingCompactor[_dBG[{{type_bundle.params}}]]] _this
    cdef public EventNotifier Notifier
{% endfor %}


{% endblock code %}
