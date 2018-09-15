{# boink/templates/processors.tpl.pyx
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
from boink.processors cimport *


cdef class FileProcessor_Base(FileProcessor):
    pass


{% for type_bundle in type_bundles %}

cdef class FileConsumer_{{type_bundle.suffix}}(FileProcessor_Base):

    def __cinit__(self, dBG_{{type_bundle.suffix}} graph,
                        uint32_t output_interval):
        self._this = make_unique[_FileConsumer[_dBG[{{type_bundle.params}}]]](graph._this.get(),
                                                                              output_interval)
        self.storage_type = graph.storage_type
        self.shifter_type = graph.shifter_type

    def process(self, str input_filename):
        deref(self._this).process(_bstring(input_filename))

        return (deref(self._this).n_reads(),
                deref(self._this).n_consumed())


cdef class DecisionNodeProcessor_{{type_bundle.suffix}}(FileProcessor_Base):
    
    def __cinit__(self, StreamingCompactor_{{type_bundle.suffix}} compactor, str output_filename,
                        uint32_t output_interval):
        self.output_filename = output_filename
        cdef string _output_filename = _bstring(output_filename)
        self._this = make_unique[_DecisionNodeProcessor[_dBG[{{type_bundle.params}}]]](compactor._this.get(),
                                                                                       _output_filename,
                                                                                       output_interval)
        self.storage_type = compactor.storage_type
        self.shifter_type = compactor.shifter_type

    def process(self, str input_filename):
        deref(self._this).process(_bstring(input_filename))

        return deref(self._this).n_reads()


cdef class StreamingCompactorProcessor_{{type_bundle.suffix}}(FileProcessor_Base):
    
    def __cinit__(self, StreamingCompactor_{{type_bundle.suffix}} compactor,
                        uint32_t output_interval):

        self._this = make_unique[_StreamingCompactorProcessor[_dBG[{{type_bundle.params}}]]](compactor._this.get(),
                                                                                             output_interval)
        self.Notifier = EventNotifier._wrap(<_EventNotifier*>self._this.get())

        self.storage_type = compactor.storage_type
        self.shifter_type = compactor.shifter_type

    def process(self, str input_filename, str right_filename=None):
        if right_filename is None:
            deref(self._this).process(_bstring(input_filename))
        else:
            deref(self._this).process(_bstring(input_filename),
                                      _bstring(right_filename))

        return deref(self._this).n_reads()
{% endfor %}


cdef object _make_file_consumer(dBG_Base graph, int output_interval):
    {% for type_bundle in type_bundles %}
    if graph.storage_type == "{{type_bundle.storage_type}}" and \
       graph.shifter_type == "{{type_bundle.shifter_type}}":
        return FileConsumer_{{type_bundle.suffix}}(graph, output_interval)
    {% endfor %}

    raise TypeError("Invalid dBG type.")


cdef object _make_decision_node_processor(StreamingCompactor_Base compactor, str filename, int output_interval):
    {% for type_bundle in type_bundles %}
    if compactor.storage_type == "{{type_bundle.storage_type}}" and \
       compactor.shifter_type == "{{type_bundle.shifter_type}}":
        return DecisionNodeProcessor_{{type_bundle.suffix}}(compactor, filename, output_interval)
    {% endfor %}

    raise TypeError("Invalid dBG type.")


cdef object _make_streaming_compactor_processor(StreamingCompactor_Base compactor, int output_interval):
    {% for type_bundle in type_bundles %}
    if compactor.storage_type == "{{type_bundle.storage_type}}" and \
       compactor.shifter_type == "{{type_bundle.shifter_type}}":
        return StreamingCompactorProcessor_{{type_bundle.suffix}}(compactor, output_interval)
    {% endfor %}

    raise TypeError("Invalid dBG type.")

{% endblock code %}
