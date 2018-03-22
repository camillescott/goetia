from libcpp.memory cimport shared_ptr, make_shared

from boink.dbg cimport *

cdef class Assembler_Base:
    cdef readonly object storage_type
    cdef readonly object shifter_type

{% for Storage_t in Storage_types %}
{% for Shifter_t in Shifter_types %}

{% set tparams %}{{Storage_t}},{{Shifter_t}}{% endset %}
{% set suffix %}{{Storage_t}}_{{Shifter_t}}{% endset %}

cdef class Assembler_{{suffix}}(Assembler_Base):
    cdef shared_ptr[_AssemblerMixin[_dBG[{{tparams}}]]] _this
    cdef shared_ptr[_dBG[{{tparams}}]] _graph
    cdef readonly dBG_{{suffix}} Graph

{% endfor %}
{% endfor %}

cdef object _make_assembler(dBG_Base graph)
