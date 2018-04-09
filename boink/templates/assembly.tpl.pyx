{# boink/templates/assembly.tpl.pyx
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}
{% extends "base.tpl" %}
{% from "dbg_types.tpl" import iter_types %}
{% block code %}

from cython.operator cimport dereference as deref

from libc.stdint cimport uint64_t
from libcpp.memory cimport make_unique
from libcpp.string cimport string

from boink.dbg cimport *
from boink.utils cimport *


cdef class Assembler_Base:
    pass

{% call(Storage_t, Shifter_t, tparams, suffix) iter_types(Storage_types, Shifter_types) %}
cdef class Assembler_{{suffix}}(Assembler_Base):

    def __cinit__(self, dBG_{{suffix}} graph):
        if type(self) is Assembler_{{suffix}}:
            self._this = make_unique[_AssemblerMixin[_dBG[{{tparams}}]]](graph._this.get())
            self._graph = graph._this.get()
            self.Graph = graph
        self.storage_type = graph.storage_type
        self.shifter_type = graph.shifter_type

    @property
    def cursor(self):
        return deref(self._this).get_cursor()

    @cursor.setter
    def cursor(self, str seed):
        deref(self._this).set_cursor(_bstring(seed))

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

    def assemble_left(self, object seed=None):
        cdef Path path
        if seed is None:
            deref(self._this).assemble_left(path)
        else:
            deref(self._this).assemble_left(_bstring(seed), path)

        return deref(self._this).to_string(path)
        
    def assemble_right(self, object seed=None):
        cdef Path path
        if seed is None:
            deref(self._this).assemble_right(path)
        else:
            deref(self._this).assemble_right(_bstring(seed), path)

        return deref(self._this).to_string(path)
{% endcall %}


cdef object _make_assembler(dBG_Base graph):
    {% call(Storage_t, Shifter_t, tparams, suffix) iter_types(Storage_types, Shifter_types) %}
    if graph.storage_type == "{{Storage_t}}" and \
       graph.shifter_type == "{{Shifter_t}}":
        return Assembler_{{suffix}}(graph)
    {% endcall %}

    raise TypeError("Invalid dBG type.")

{% endblock code %}
