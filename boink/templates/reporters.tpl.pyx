{# boink/templates/reporters.tpl.pyx
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}
{% extends "base.tpl" %}
{% block code %}

from boink.dbg cimport *
from boink.compactor cimport *
from boink.utils cimport *

import sys

from cython.operator cimport dereference as deref
from libcpp.memory cimport unique_ptr, make_unique

from boink.utils cimport _bstring
from boink.events cimport _EventListener


cdef class StreamingCompactorReporter_Base(SingleFileReporter):
    pass

{% for type_bundle in type_bundles %}

cdef class StreamingCompactorReporter_{{type_bundle.suffix}}(StreamingCompactorReporter_Base):
    
    def __cinit__(self, str output_filename, StreamingCompactor_{{type_bundle.suffix}} compactor,
                        *args, **kwargs):
        if type(self) is StreamingCompactorReporter_{{type_bundle.suffix}}:
            self._s_owner = make_unique[_StreamingCompactorReporter[_dBG[{{type_bundle.params}}]]](\
                    compactor._this.get(), _bstring(output_filename))
            self._s_this = self._s_owner.get()
            self._this = self._s_this
            self._listener = <_EventListener*>self._s_owner.get()

        self.storage_type = compactor.storage_type
        self.shifter_type = compactor.shifter_type

cdef class cDBGWriter_{{type_bundle.suffix}}(cDBGWriter_Base):

    def __cinit__(self, str output_prefix, str graph_format, cDBG_{{type_bundle.suffix}} cdbg,
                        *args, **kwargs):
        
        self.storage_type = cdbg.storage_type
        self.shifter_type = cdbg.shifter_type

        if type(self) is cDBGWriter_{{type_bundle.suffix}}:
            self._s_owner = make_unique[_cDBGWriter[_dBG[{{type_bundle.params}}]]](cdbg._this,
                                                     convert_format(graph_format),
                                                     _bstring(output_prefix))
            self._s_this = self._s_owner.get()
            self._this = self._s_this
            self._listener = <_EventListener*>self._s_owner.get()

{% endfor %}

cdef object _make_streaming_compactor_reporter(str output_filename,
                                               StreamingCompactor_Base compactor):
    {% for type_bundle in type_bundles %}
    if compactor.storage_type == "{{type_bundle.storage_type}}" and \
       compactor.shifter_type == "{{type_bundle.shifter_type}}":
        return StreamingCompactorReporter_{{type_bundle.suffix}}(output_filename, compactor)
    {% endfor %}

    raise TypeError("Invalid dBG type.")

cdef object _make_cdbgwriter_reporter(str output_prefix,
                                      str graph_format,
                                      cDBG_Base cdbg):
    {% for type_bundle in type_bundles %}
    if cdbg.storage_type == "{{type_bundle.storage_type}}" and \
       cdbg.shifter_type == "{{type_bundle.shifter_type}}":
        return cDBGWriter_{{type_bundle.suffix}}(output_prefix, graph_format, cdbg)
    {% endfor %}

    raise TypeError("Invalid dBG type.")
{% endblock code %}
