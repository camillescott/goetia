{# boink/templates/assembly.tpl.pyx
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}
{% extends "base.tpl" %}
{% block code %}

from cython.operator cimport dereference as deref

from libc.stdint cimport uint64_t
from libcpp.memory cimport make_shared
from libcpp.string cimport string

from boink.dbg cimport *
from boink.utils cimport *


cdef class Assembler:

    @staticmethod
    def build(dBG graph):
        {% for type_bundle in type_bundles %}
        if graph.storage_type == "{{type_bundle.storage_type}}" and \
           graph.shifter_type == "{{type_bundle.shifter_type}}":
            return Assembler_{{type_bundle.suffix}}(graph)
        {% endfor %}

        raise TypeError("Invalid dBG type.")


{% for type_bundle in type_bundles %}
cdef class Assembler_{{type_bundle.suffix}}(Assembler):

    def __cinit__(self, dBG_{{type_bundle.suffix}} graph):
        if type(self) is Assembler_{{type_bundle.suffix}}:
            self._this = make_shared[_AssemblerMixin[_dBG[{{type_bundle.params}}]]](graph._this)
            self._graph = graph._this
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
{% endfor %}

{% endblock code %}
