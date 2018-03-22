from libcpp.memory cimport shared_ptr, make_shared

cdef class dBG_Base:
    cdef readonly object storage_type
    cdef readonly object shifter_type

{% for Storage_t in Storage_types %}
{% for Shifter_t in Shifter_types %}
{% set tparams %}{{Storage_t}},{{Shifter_t}}{% endset %}
{% set suffix %}{{Storage_t}}_{{Shifter_t}}{% endset %}

cdef class dBG_{{suffix}}(dBG_Base):
    cdef shared_ptr[_dBG[{{tparams}}]] _this
    cdef hash_t _handle_kmer(self, object) except 0
{% endfor %}
{% endfor %}

cdef object _make_dbg(int K, uint64_t starting_size, int n_tables, 
                     str storage=*, str shifter=*)
