{# boink/templates/cdbg.tpl.pxd
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}
{% extends "base.tpl" %}
{% block code %}


cdef class cDBG_Base:
    cdef readonly object shifter_type
    cdef readonly object storage_type

{% for type_bundle in type_bundles %}
cdef class cDBG_{{type_bundle.suffix}}(cDBG_Base):
    cdef _cDBG[_dBG[{{type_bundle.params}}]] * _this
    cdef public EventNotifier Notifier

    @staticmethod
    cdef cDBG_{{type_bundle.suffix}} _wrap(_cDBG[_dBG[{{type_bundle.params}}]] *)
{% endfor %}

{% endblock code %}
