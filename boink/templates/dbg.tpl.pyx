{# boink/templates/dbg.tpl.pyx
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
from libcpp.vector cimport vector

from khmer._oxli.utils import get_n_primes_near_x, is_str, is_num
from boink.utils cimport _bstring, _ustring


cdef class dBG_Base:
    
    def __cinit__(self, *args, **kwargs):
        self.allocated = False

{% call(Storage_t, Shifter_t, tparams, suffix) iter_types(Storage_types, Shifter_types) %}
cdef class dBG_{{suffix}}(dBG_Base):

    def __cinit__(self, int K, uint64_t starting_size, int n_tables,
                  *args, **kwargs):
        cdef vector[uint64_t] primes
        #if type(self) is dBG_{{suffix}}:
        if not self._this:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self._this = make_unique[_dBG[{{tparams}}]](K, primes)
            self._assembler = make_unique[_AssemblerMixin[_dBG[{{tparams}}]]](self._this.get())
            self.allocated = True

        self.storage_type = "{{Storage_t}}"
        self.shifter_type = "{{Shifter_t}}"
        self.suffix = "{{suffix}}"

    cdef hash_t _handle_kmer(self, object kmer) except 0:
        cdef hash_t handled
        if is_num(kmer):
            handled = <hash_t> kmer
        else:
            handled = deref(self._this).hash(_bstring(kmer))
        return handled

    def clone(self):
        cdef dBG_{{suffix}} cloned = dBG_{{suffix}}(1,1,1)
        cloned._this = deref(self._this).clone()
        cloned._assembler.reset(new _AssemblerMixin[_dBG[{{tparams}}]](self._this.get()))
        return cloned

    def add(self, object kmer):
        return deref(self._this).add(self._handle_kmer(kmer))

    def get(self, object kmer):
        return deref(self._this).get(self._handle_kmer(kmer))

    def hash(self, str kmer):
        return deref(self._this).hash(_bstring(kmer))

    def hashes(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef unique_ptr[_KmerIterator[{{Shifter_t}}]] kmer_iter = \
                deref(self._this).get_hash_iter(_sequence)

        cdef hash_t h
        while(not deref(kmer_iter).done()):
            h = deref(kmer_iter).next()
            yield h

    def add_sequence_and_report(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef vector[bool] report
        deref(self._this).add_sequence(_sequence, report)
        return report

    def add_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        return deref(self._this).add_sequence(_sequence)

    # compatibility with oxli API
    consume = add_sequence

    def get_counts(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef list counts = deref(self._this).get_counts(_sequence)
        return counts

    def left_degree(self, str kmer):
        deref(self._assembler).set_cursor(_bstring(kmer))
        return deref(self._assembler).degree_left()

    def right_degree(self, str kmer):
        deref(self._assembler).set_cursor(_bstring(kmer))
        return deref(self._assembler).degree_right()

    def degree(self, str kmer):
        deref(self._assembler).set_cursor(_bstring(kmer))
        return deref(self._assembler).degree_left() + \
               deref(self._assembler).degree_right()

    @property
    def n_unique(self):
        return deref(self._this).n_unique()

    @property
    def n_occupied(self):
        return deref(self._this).n_occupied()

    @property
    def K(self):
        return deref(self._this).K()

    def save(self, file_name):
        deref(self._this).save(_bstring(file_name))

    @classmethod
    def load(cls, file_name):
        cdef dBG_{{suffix}} obj = cls(1, 1, 1)
        deref(obj._this).load(_bstring(file_name))
        return obj
{% endcall %}

cdef object _make_dbg(int K, uint64_t starting_size, int n_tables,
                      str storage='BitStorage',
                      str shifter='DefaultShifter'):

    {% call(Storage_t, Shifter_t, tparams, suffix) iter_types(Storage_types, Shifter_types) %}
    if storage == "{{Storage_t}}" and shifter == "{{Shifter_t}}":
        return dBG_{{suffix}}(K, starting_size, n_tables)
    {% endcall %}
    raise TypeError("Invalid Storage or Shifter type.")


def get_dbg_type(str storage='BitStorage',
                 str shifter='DefaultShifter'):
    {% call(Storage_t, Shifter_t, tparams, suffix) iter_types(Storage_types, Shifter_types) %}
    if storage == "{{Storage_t}}" and shifter == "{{Shifter_t}}":
        return dBG_{{suffix}}
    {% endcall %}
    raise TypeError("Invalid Storage or Shifter type: ({0},{1})".format(storage, shifter))


{% endblock code %}
