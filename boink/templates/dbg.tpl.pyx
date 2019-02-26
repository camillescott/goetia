{# boink/templates/dbg.tpl.pyx
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
from libcpp.vector cimport vector

from boink.utils cimport (_bstring, _ustring, get_n_primes_near_x,
                          is_str, is_num)


cdef class dBG:
    
    def __cinit__(self, *args, **kwargs):
        self.allocated = False

    @staticmethod
    def build(int K,
              *args,
              #uint64_t starting_size,
              #int n_tables,
              str storage='_BitStorage',
              str shifter='_DefaultShifter',
              **kwargs):
        {% for type_bundle in type_bundles %}
        if storage == "{{type_bundle.storage_type}}" and \
           shifter == "{{type_bundle.shifter_type}}":
            return dBG_{{type_bundle.suffix}}(K, *args, **kwargs)

        {% endfor %}
        raise TypeError("Invalid Storage or Shifter type: ({0},{1})".format(storage, shifter))

    @staticmethod
    def get_type(str storage='_BitStorage',
                 str shifter='_DefaultShifter'):
        {% for type_bundle in type_bundles %}
        if storage == "{{type_bundle.storage_type}}" and \
           shifter == "{{type_bundle.shifter_type}}":
            return dBG_{{type_bundle.suffix}}
        {% endfor %}
        raise TypeError("Invalid Storage or Shifter type: ({0},{1})".format(storage, shifter))


{% for type_bundle in type_bundles %}
cdef class dBG_{{type_bundle.suffix}}(dBG):

    def __cinit__(self, int K,# uint64_t starting_size, int n_tables,
                  *args, **kwargs):
        #if type(self) is dBG_{{suffix}}:
        self.storage_type = "{{type_bundle.storage_type}}"
        self.shifter_type = "{{type_bundle.shifter_type}}"
        self.suffix = "{{type_bundle.suffix}}"

        cdef int starting_size
        cdef int n_tables
        if self.storage_type in ['_BitStorage', '_ByteStorage', '_NibbleStorage']:
            starting_size, n_tables = args

        if not self._this:
            {% if type_bundle.storage_type  in ['_BitStorage', '_ByteStorage', '_NibbleStorage'] %}
            self._this = make_shared[_dBG[{{type_bundle.params}}]](K, starting_size, n_tables)
            {% else %}
            self._this = make_shared[_dBG[{{type_bundle.params}}]](K)
            {% endif %}

            self._assembler = make_shared[_AssemblerMixin[_dBG[{{type_bundle.params}}]]](self._this)
            self.allocated = True

    cdef hash_t _handle_kmer(self, object kmer) except 0:
        cdef hash_t handled
        if is_num(kmer):
            handled = <hash_t> kmer
        else:
            handled = deref(self._this).hash(_bstring(kmer))
        return handled

    def insert(self, object kmer):
        return deref(self._this).insert(self._handle_kmer(kmer))

    def insert_and_query(self, object kmer):
        return deref(self._this).insert_and_query(self._handle_kmer(kmer))


    def query(self, object kmer):
        return deref(self._this).query(self._handle_kmer(kmer))

    def hash(self, str kmer):
        return deref(self._this).hash(_bstring(kmer))

    def hashes(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef shared_ptr[_KmerIterator[{{type_bundle.shifter_type}}]] kmer_iter = \
                deref(self._this).get_hash_iter(_sequence)

        cdef hash_t h
        while(not deref(kmer_iter).done()):
            h = deref(kmer_iter).next()
            yield h

    def neighbors(self, str root):
        cdef bytes _root = _bstring(root)
        cdef pair[vector[kmer_t], vector[kmer_t]] result = deref(self._this).neighbor_kmers(_root)
        
        cdef list left = []
        cdef list right = []
        for i in range(result.first.size()):
            left.append((result.first[i].hash, result.first[i].kmer))
        for i in range(result.second.size()):
            right.append((result.second[i].hash, result.second[i].kmer))

        return left, right

    def insert_sequence_and_report(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef vector[hash_t] _hashes
        cdef vector[count_t] _report
        deref(self._this).insert_sequence(_sequence,
                                       _hashes,
                                       _report)
        return _hashes, _report

    def insert_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        return deref(self._this).insert_sequence(_sequence)

    # compatibility with oxli API
    consume = insert_sequence

    def query_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef list counts = deref(self._this).query_sequence(_sequence)
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
        cdef dBG_{{type_bundle.suffix}} obj = cls(1, 1, 1)
        deref(obj._this).load(_bstring(file_name))
        return obj

    def clone(self):
        cdef dBG_{{type_bundle.suffix}} cloned = dBG_{{type_bundle.suffix}}(1,1,1)
        cloned._this = deref(self._this).clone()
        cloned._assembler = make_shared[_AssemblerMixin[_dBG[{{type_bundle.params}}]]](self._this)
        return cloned

    def reset(self):
        deref(self._this).reset()
{% endfor %}

{% endblock code %}
