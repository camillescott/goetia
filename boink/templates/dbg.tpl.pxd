{# boink/templates/dbg.tpl.pxd
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}
{% extends "base.tpl" %}
{% block code %}

from libcpp.memory cimport unique_ptr
from boink.assembly cimport _AssemblerMixin

cdef class dBG_Base:
    cdef readonly object storage_type
    cdef readonly object shifter_type
    cdef readonly object suffix
    cdef object allocated

{% for type_bundle in type_bundles %}
cdef class dBG_{{type_bundle.suffix}}(dBG_Base):
    cdef unique_ptr[_dBG[{{type_bundle.params}}]] _this
    cdef unique_ptr[_AssemblerMixin[_dBG[{{type_bundle.params}}]]] _assembler
    cdef hash_t _handle_kmer(self, object) except 0
{% endfor %}

cdef object _make_dbg(int K, uint64_t starting_size, int n_tables, 
                      str storage=*, str shifter=*)

{% endblock code %}
