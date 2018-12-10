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


cdef class StreamingCompactor_Base:
    cdef readonly object storage_type
    cdef readonly object shifter_type
    cdef object graph

{% for type_bundle in type_bundles %}
cdef class StreamingCompactor_{{type_bundle.suffix}}(StreamingCompactor_Base):
    cdef shared_ptr[_StreamingCompactor[_dBG[{{type_bundle.params}}]]] _this
    cdef shared_ptr[_dBG[{{type_bundle.params}}]] _graph
    cdef public cDBG_{{type_bundle.suffix}} cdbg
    cdef public EventNotifier Notifier
{% endfor %}

cdef object _make_streaming_compactor(dBG_Base graph, Instrumentation inst=*)

{% endblock code %}
