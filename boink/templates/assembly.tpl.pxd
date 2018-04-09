{# boink/templates/assembly.tpl.pxd
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}
{% extends "base.tpl" %}
{% from "dbg_types.tpl" import iter_types %}
{% block code %}

from libcpp.memory cimport unique_ptr

from boink.dbg cimport *


cdef class Assembler_Base:
    cdef readonly object storage_type
    cdef readonly object shifter_type


{% call(Storage_t, Shifter_t, tparams, suffix) iter_types(Storage_types, Shifter_types) %}
cdef class Assembler_{{suffix}}(Assembler_Base):
    cdef unique_ptr[_AssemblerMixin[_dBG[{{tparams}}]]] _this
    cdef _dBG[{{tparams}}] * _graph
    cdef readonly dBG_{{suffix}} Graph
{% endcall %}

cdef object _make_assembler(dBG_Base graph)

{% endblock code %}
