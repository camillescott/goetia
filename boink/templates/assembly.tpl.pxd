{# boink/templates/assembly.tpl.pxd
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


cdef class Assembler_Base:
    cdef readonly object storage_type
    cdef readonly object shifter_type

{% for type_bundle in type_bundles %}
cdef class Assembler_{{type_bundle.suffix}}(Assembler_Base):
    cdef shared_ptr[_AssemblerMixin[_dBG[{{type_bundle.params}}]]] _this
    cdef shared_ptr[_dBG[{{type_bundle.params}}]] _graph
    cdef readonly dBG_{{type_bundle.suffix}} Graph
{% endfor %}

cdef object _make_assembler(dBG graph)

{% endblock code %}
