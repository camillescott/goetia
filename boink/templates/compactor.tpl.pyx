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
from libcpp.memory cimport make_unique
from libcpp.string cimport string

from boink.dbg cimport *
from boink.utils cimport *


cdef class StreamingCompactor_Base:
    pass

{% for type_bundle in type_bundles %}
cdef class StreamingCompactor_{{type_bundle.suffix}}(StreamingCompactor_Base):

    def __cinit__(self, dBG_{{type_bundle.suffix}} graph):

        if type(self) is StreamingCompactor_{{type_bundle.suffix}}:
            self._this = make_unique[_StreamingCompactor[_dBG[{{type_bundle.params}}]]](graph._this.get())
            self._graph = graph._this.get()
            self.cdbg = cDBG._wrap(deref(self._this).cdbg)
            self.Notifier = EventNotifier._wrap(<_EventNotifier*>self._this.get())
        
        self.storage_type = graph.storage_type
        self.shifter_type = graph.shifter_type

    def find_decision_kmers(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        cdef vector[uint32_t] positions
        cdef vector[hash_t] hashes
        cdef vector[NeighborBundle] neighbors

        deref(self._this).find_decision_kmers(_sequence,
                                               positions,
                                               hashes,
                                               neighbors)

        return positions, hashes

    def update_sequence(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        deref(self._this).update_sequence(_sequence)

    def find_new_segments(self, str sequence):
        cdef string _sequence = _bstring(sequence)

        cdef deque[_compact_segment] _segments
        deref(self._this).find_new_segments(_sequence,
                                            _segments)

        segments = []
        cdef size_t i = 0
        for i in range(_segments.size()):
            if _segments[i].is_null():
                segment = Segment(sequence = '',
                                  is_decision_kmer = False,
                                  left_anchor = 0,
                                  right_anchor = 0,
                                  start = 0,
                                  length = 0,
                                  is_null = True)
            else:
                segment_seq = sequence[_segments[i].start_pos : \
                                       _segments[i].start_pos + _segments[i].length]
                segment = Segment(sequence =           segment_seq,
                                  is_decision_kmer =   _segments[i].is_decision_kmer,
                                  left_anchor =        _segments[i].left_anchor,
                                  right_anchor =       _segments[i].right_anchor,
                                  start =              _segments[i].start_pos,
                                  length =             _segments[i].length,
                                  is_null =            False)
            segments.append(segment)

        return segments

{% endfor %}


cdef object _make_streaming_compactor(dBG_Base graph):
    {% for type_bundle in type_bundles %}
    if graph.storage_type == "{{type_bundle.storage_type}}" and \
       graph.shifter_type == "{{type_bundle.shifter_type}}":
        return StreamingCompactor_{{type_bundle.suffix}}(graph)
    {% endfor %}

    raise TypeError("Invalid dBG type.")

{% endblock code %}
