{# boink/templates/compactor.tpl.pxd
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}
{% extends "base.tpl" %}
{% block code %}

from libcpp.memory cimport unique_ptr

from boink.dbg cimport *
from boink.cdbg cimport *


cdef class StreamingCompactor_Base:
    cdef readonly object storage_type
    cdef readonly object shifter_type

{% for type_bundle in type_bundles %}
cdef class StreamingCompactor_{{type_bundle.suffix}}(StreamingCompactor_Base):
    cdef unique_ptr[_StreamingCompactor[_dBG[{{type_bundle.params}}]]] _this
    cdef _dBG[{{type_bundle.params}}] * _graph
    cdef public cDBG cdbg
    cdef public EventNotifier Notifier
{% endfor %}

cdef object _make_streaming_compactor(dBG_Base graph)

{% endblock code %}
