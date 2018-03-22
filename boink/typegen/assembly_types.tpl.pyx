from cython.operator cimport dereference as deref

from libc.stdint cimport uint64_t
from libcpp.string cimport string

from boink import dbg
from boink.utils cimport *

cdef class Assembler_Base:
    pass


{% for Storage_t in Storage_types %}
{% for Shifter_t in Shifter_types %}
{% set tparams %}{{Storage_t}},{{Shifter_t}}{% endset %}
{% set suffix %}{{Storage_t}}_{{Shifter_t}}{% endset %}

cdef class Assembler_{{suffix}}(Assembler_Base):

    def __cinit__(self, dBG_{{suffix}} graph):
        if type(self) is Assembler_{{suffix}}:
            self._this = make_shared[_AssemblerMixin[_dBG[{{tparams}}]]](graph._this)
            self._graph = graph._this
            self.Graph = graph
        self.storage_type = graph.storage_type
        self.shifter_type = graph.shifter_type

    @property
    def cursor(self):
        return deref(self._this).get_cursor()

    def clear_seen(self):
        deref(self._this).clear_seen()

    def degree_left(self):
        return deref(self._this).degree_left()

    def degree_right(self):
        return deref(self._this).degree_right()

    def degree(self):
        return deref(self._this).degree()

    def assemble(self, str seed):
        cdef bytes _seed = _bstring(seed)
        cdef Path path

        deref(self._this).assemble(_seed, path)
        return deref(self._this).to_string(path)

{% endfor %}
{% endfor %}

cdef object _make_assembler(dBG_Base graph):
    {% for Storage_t in Storage_types %}
    {% set outer_first = loop.first %}
    {% for Shifter_t in Shifter_types %}
    {% set conditional %}graph.storage_type == "{{Storage_t}}" and graph.shifter_type == "{{Shifter_t}}"{% endset %}
    {% set suffix %}{{Storage_t}}_{{Shifter_t}}{% endset %}
    {% if outer_first and loop.first %}
    if {{conditional}}:
    {% else %}
    elif {{conditional}}:
    {% endif %}
        return Assembler_{{suffix}}(graph)
    {% endfor %}
    {% endfor %}
    else:
        raise TypeError("Invalid dBG type.")
