{# boink/templates/dbg.tpl.pxd
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}
{% extends "base.tpl" %}
{% from "dbg_types.tpl" import iter_types %}
{% block code %}

from libcpp.memory cimport shared_ptr, make_shared
from boink.assembly cimport _AssemblerMixin

cdef class dBG_Base:
    cdef readonly object storage_type
    cdef readonly object shifter_type
    cdef object allocated


{% call(Storage_t, Shifter_t, tparams, suffix) iter_types(Storage_types, Shifter_types) %}
cdef class dBG_{{suffix}}(dBG_Base):
    cdef shared_ptr[_dBG[{{tparams}}]] _this
    cdef shared_ptr[_AssemblerMixin[_dBG[{{tparams}}]]] _assembler
    cdef hash_t _handle_kmer(self, object) except 0
{% endcall %}

cdef object _make_dbg(int K, uint64_t starting_size, int n_tables, 
                     str storage=*, str shifter=*)

{% endblock code %}
