{# boink/templates/cdbg.tpl.pyx
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}
{% extends "base.tpl" %}
{% block code %}

cdef class cDBG_Base:

    SAVE_FORMATS = ['graphml', 'edgelist',
                    'gfa1', 'gfa2', 'fasta', 'gml']

    def __cinit__(self, *args, **kwargs):
        pass


cdef cDBGFormat convert_format(str file_format) except *:
    if file_format in cDBG_Base.SAVE_FORMATS:
        if file_format == 'graphml':
            return cDBGFormat.GRAPHML
        elif file_format == 'fasta':
            return cDBGFormat.FASTA
        elif file_format == 'gfa1':
            return cDBGFormat.GFA1
        else:
            raise NotImplementedError("Support for {0} not yet "
                                      "implemented".format(file_format))
    else:
        formats = ', '.join(cDBG_Base.SAVE_FORMATS)
        raise ValueError("{0} not a valid save format. "
                         "Format must be one of: {1}".format(file_format,
                                                             formats))

{% for type_bundle in type_bundles %}
cdef class cDBG_{{type_bundle.suffix}}(cDBG_Base):

    @staticmethod
    cdef cDBG_{{type_bundle.suffix}} _wrap(_cDBG[_dBG[{{type_bundle.params}}]] * ptr):
        cdef cDBG_{{type_bundle.suffix}} cdbg = cDBG_{{type_bundle.suffix}}()
        cdbg._this = ptr
        cdbg.Notifier = EventNotifier._wrap(<_EventNotifier*>cdbg._this)
        cdbg.storage_type = "{{type_bundle.storage_type}}"
        cdbg.shifter_type = "{{type_bundle.shifter_type}}"
        return cdbg
    '''
    def unodes(self):
        cdef _cDBG[_dBG[{{type_bundle.params}}]].unode_iter_t it = \
                deref(self._this).unodes_begin()
        cdef _UnitigNode * unode
        while(it != deref(self._this).unodes_end()):
            unode = deref(it).second.get()
            yield UnitigNodeView._wrap(unode)
            princ(it)


    def dnodes(self):
        cdef _cDBG[_dBG[{{type_bundle.params}}]].dnode_iter_t it = \
                deref(self._this).dnodes_begin()
        cdef _DecisionNode * dnode
        while(it != deref(self._this).dnodes_end()):
            dnode = deref(it).second.get()
            yield DecisionNodeView._wrap(dnode)
            princ(it)
    '''

    def query_unode_hash(self, hash_t h):
        cdef _UnitigNode * _node = deref(self._this).query_unode_tag(h)
        if _node != NULL:
            return UnitigNodeView._wrap(_node)
        else:
            return None

    def query_unode_id(self, id_t node_id):
        cdef _UnitigNode * _node = \
            deref(self._this).query_unode_id(node_id)
        if _node != NULL:
            return UnitigNodeView._wrap(_node)
        else:
            return None

    def query_unode_end(self, hash_t h):
        cdef _UnitigNode * _node = \
            deref(self._this).query_unode_end(h)
        if _node != NULL:
            return UnitigNodeView._wrap(_node)
        else:
            return None
    
    def query_dnode(self, hash_t h):
        cdef _DecisionNode * _node = \
            deref(self._this).query_dnode(h)
        if _node != NULL:
            return DecisionNodeView._wrap(_node)
        else:
            return None

    def has_dnode(self, hash_t h):
        return deref(self._this).has_dnode(h)

    def query_dnodes(self, str sequence):
        pass

    @property
    def n_updates(self):
        return deref(self._this).n_updates()

    @property
    def n_unodes(self):
        return deref(self._this).n_unitig_nodes()

    @property
    def n_dnodes(self):
        return deref(self._this).n_decision_nodes()

    @property
    def n_tags(self):
        return deref(self._this).n_tags()

    @property
    def n_unitig_ends(self):
        return deref(self._this).n_unitig_ends()

    def validate(self, str filename):
        deref(self._this).validate(_bstring(filename))

    def save(self, str filename, str file_format):
        if file_format is None:
            return
        else:
            deref(self._this).write(_bstring(filename),
                                    convert_format(file_format))

{% endfor %}

def get_cdbg_type(str storage='_BitStorage',
                  str shifter='_DefaultShifter'):
    {% for type_bundle in type_bundles %}
    if storage == "{{type_bundle.storage_type}}" and \
       shifter == "{{type_bundle.shifter_type}}":
        return cDBG_{{type_bundle.suffix}}
    {% endfor %}
    raise TypeError("Invalid Storage or Shifter type: ({0},{1})".format(storage, shifter))


{% endblock code %}
