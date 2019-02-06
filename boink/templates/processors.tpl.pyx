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
from libcpp.memory cimport make_shared
from libcpp.string cimport string

from boink.dbg cimport *
from boink.utils cimport *
from boink.processors cimport *


cdef class FileConsumer(FileProcessor):

    @staticmethod
    def build(dBG graph,
              uint64_t fine_interval,
              uint64_t medium_interval,
              uint64_t coarse_interval):
    
        {% for type_bundle in type_bundles %}
        if graph.storage_type == "{{type_bundle.storage_type}}" and \
           graph.shifter_type == "{{type_bundle.shifter_type}}":
            return FileConsumer_{{type_bundle.suffix}}(graph, 
                                                       fine_interval,
                                                       medium_interval,
                                                       coarse_interval)
        {% endfor %}

        raise TypeError("Invalid dBG type.")


cdef class DecisionNodeProcessor(FileProcessor):
    
    @staticmethod
    def build(StreamingCompactor compactor,
              str filename, 
              uint64_t fine_interval,
              uint64_t medium_interval,
              uint64_t coarse_interval):

        {% for type_bundle in type_bundles %}
        if compactor.storage_type == "{{type_bundle.storage_type}}" and \
           compactor.shifter_type == "{{type_bundle.shifter_type}}":
            return DecisionNodeProcessor_{{type_bundle.suffix}}(compactor,
                                                                filename,  
                                                                fine_interval,
                                                                medium_interval,
                                                                coarse_interval)
        {% endfor %}

        raise TypeError("Invalid dBG type.")


cdef class StreamingCompactorProcessor(FileProcessor):

    @staticmethod
    def build(StreamingCompactor compactor, 
              uint64_t fine_interval,
              uint64_t medium_interval,
              uint64_t coarse_interval):

        {% for type_bundle in type_bundles %}
        if compactor.storage_type == "{{type_bundle.storage_type}}" and \
           compactor.shifter_type == "{{type_bundle.shifter_type}}":
            return StreamingCompactorProcessor_{{type_bundle.suffix}}(compactor,
                                                                      fine_interval,
                                                                      medium_interval,
                                                                      coarse_interval)
        {% endfor %}

        raise TypeError("Invalid dBG type.")


cdef class NormalizingCompactor(FileProcessor):

    @staticmethod
    def build(StreamingCompactor compactor, 
              unsigned int cutoff,
              uint64_t fine_interval,
              uint64_t medium_interval,
              uint64_t coarse_interval):

        {% for type_bundle in type_bundles %}
        if compactor.storage_type == "{{type_bundle.storage_type}}" and \
           compactor.shifter_type == "{{type_bundle.shifter_type}}":
            return NormalizingCompactor_{{type_bundle.suffix}}(compactor,
                                                               cutoff,
                                                               fine_interval,
                                                               medium_interval,
                                                               coarse_interval)
        {% endfor %}

        raise TypeError("Invalid dBG type.")


{% for type_bundle in type_bundles %}

cdef class FileConsumer_{{type_bundle.suffix}}(FileConsumer):

    def __cinit__(self, dBG_{{type_bundle.suffix}} graph,
                        uint64_t fine_interval,
                        uint64_t medium_interval,
                        uint64_t coarse_interval):

        self._this = make_shared[_FileConsumer[_dBG[{{type_bundle.params}}]]](graph._this,
                                                                              fine_interval,
                                                                              medium_interval,
                                                                              coarse_interval)
        self.storage_type = graph.storage_type
        self.shifter_type = graph.shifter_type

    def process(self, str input_filename):
        deref(self._this).process(_bstring(input_filename))

        return (deref(self._this).n_reads(),
                deref(self._this).n_consumed())


cdef class DecisionNodeProcessor_{{type_bundle.suffix}}(DecisionNodeProcessor):
    
    def __cinit__(self, StreamingCompactor_{{type_bundle.suffix}} compactor,
                        str output_filename,
                        uint64_t fine_interval,
                        uint64_t medium_interval,
                        uint64_t coarse_interval):

        self.output_filename = output_filename
        cdef string _output_filename = _bstring(output_filename)
        self._this = make_shared[_DecisionNodeProcessor[_dBG[{{type_bundle.params}}]]](compactor._this,
                                                                                       _output_filename,
                                                                                       fine_interval,
                                                                                       medium_interval,
                                                                                       coarse_interval)
        self.storage_type = compactor.storage_type
        self.shifter_type = compactor.shifter_type

    def process(self, str input_filename):
        deref(self._this).process(_bstring(input_filename))

        return deref(self._this).n_reads()


cdef class StreamingCompactorProcessor_{{type_bundle.suffix}}(StreamingCompactorProcessor):
    
    def __cinit__(self, StreamingCompactor_{{type_bundle.suffix}} compactor,
                        uint64_t fine_interval,
                        uint64_t medium_interval,
                        uint64_t coarse_interval):

        self._this = make_shared[_StreamingCompactorProcessor[_dBG[{{type_bundle.params}}]]](compactor._this,
                                                                                             fine_interval,
                                                                                             medium_interval,
                                                                                             coarse_interval)
        self.Notifier = EventNotifier._wrap(<shared_ptr[_EventNotifier]>self._this)

        self.storage_type = compactor.storage_type
        self.shifter_type = compactor.shifter_type

    def process(self, str input_filename, str right_filename=None):
        if right_filename is None:
            deref(self._this).process(_bstring(input_filename))
        else:
            deref(self._this).process(_bstring(input_filename),
                                      _bstring(right_filename))

        return deref(self._this).n_reads()


cdef class NormalizingCompactor_{{type_bundle.suffix}}(NormalizingCompactor):
    
    def __cinit__(self, StreamingCompactor_{{type_bundle.suffix}} compactor,
                        unsigned int cutoff,
                        uint64_t fine_interval,
                        uint64_t medium_interval,
                        uint64_t coarse_interval):

        self._this = make_shared[_NormalizingCompactor[_dBG[{{type_bundle.params}}]]](compactor._this,
                                                                                      cutoff,
                                                                                      fine_interval,
                                                                                      medium_interval,
                                                                                      coarse_interval)
        self.Notifier = EventNotifier._wrap(<shared_ptr[_EventNotifier]>self._this)

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

{% endblock code %}
