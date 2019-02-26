{# boink/templates/dbg.tpl.pxd
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}
{% extends "base.tpl" %}
{% block code %}

from libcpp.memory cimport shared_ptr
from boink.assembly cimport _AssemblerMixin

cdef class dBG:
    cdef readonly object storage_type
    cdef readonly object shifter_type
    cdef readonly object suffix
    cdef object allocated


{% for type_bundle in type_bundles %}
cdef class dBG_{{type_bundle.suffix}}(dBG):
    cdef shared_ptr[_dBG[{{type_bundle.params}}]] _this
    cdef shared_ptr[_AssemblerMixin[_dBG[{{type_bundle.params}}]]] _assembler
    cdef hash_t _handle_kmer(self, object) except 0
{% endfor %}

{% endblock code %}
