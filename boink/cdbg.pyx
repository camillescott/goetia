# cdbg.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as princ
from libcpp.memory cimport make_unique

from boink.utils cimport _bstring, _ustring, make_pair
from khmer._oxli.sequence cimport Alphabets


cdef class CompactNode:

    @staticmethod
    cdef CompactNode _wrap(_CompactNode * ptr):
        cdef CompactNode node = CompactNode()
        node._cn_this = ptr
        return node

    def _check_ptr(self):
        if (self._cn_this == NULL):
            raise ValueError("_CompactNode * is stale")

    @property
    def sequence(self):
        self._check_ptr()
        return deref(self._cn_this).sequence

    def __len__(self):
        self._check_ptr()
        return deref(self._cn_this).sequence.length()

    def __str__(self):
        self._check_ptr()
        return 'CompactNode: L={0} sequence={1}'.format(len(self), self.sequence)

    def __repr__(self):
        self._check_ptr()
        return str(self)

    @property
    def node_id(self):
        self._check_ptr()
        return deref(self._cn_this).node_id

    def __eq__(self, other):
        return self.node_id == other.node_id


cdef class DecisionNode(CompactNode):

    @staticmethod
    cdef DecisionNode _wrap(_DecisionNode * ptr):
        cdef DecisionNode node = DecisionNode()
        node._dn_this = ptr
        node._cn_this = <_CompactNode*>ptr
        return node

    def _check_ptr(self):
        if (self._dn_this == NULL):
            raise ValueError("_DecisionNode * is stale")

    @property
    def count(self):
        self._check_ptr()
        return deref(self._dn_this).count()

    @property
    def right_degree(self):
        self._check_ptr()
        return deref(self._dn_this).right_degree()

    @property
    def left_degree(self):
        self._check_ptr()
        return deref(self._dn_this).left_degree()

    @property
    def degree(self):
        self._check_ptr()
        return deref(self._dn_this).degree()

    def __str__(self):
        self._check_ptr()
        return deref(self._dn_this).repr()


cdef class UnitigNode(CompactNode):

    @staticmethod
    cdef UnitigNode _wrap(_UnitigNode * ptr):
        cdef UnitigNode node = UnitigNode()
        node._un_this = ptr
        node._cn_this = <_CompactNode*>ptr
        return node

    def _check_ptr(self):
        if (self._un_this == NULL):
            raise ValueError("_UnitigNode * is stale")

    @property
    def left_end(self):
        self._check_ptr()
        return deref(self._un_this).left_end()

    @property
    def right_end(self):
        self._check_ptr()
        return deref(self._un_this).right_end()

    @property
    def meta(self):
        return _ustring(node_meta_repr(deref(self._un_this).meta()))

    def tags(self):
        self._check_ptr()
        cdef hash_t tag
        for tag in deref(self._un_this).tags:
            yield tag

    def __str__(self):
        self._check_ptr()
        return deref(self._un_this).repr()


cdef cDBGFormat convert_format(str file_format) except *:
    if file_format in cDBG.SAVE_FORMATS:
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
        formats = ', '.join(cDBG.SAVE_FORMATS)
        raise ValueError("{0} not a valid save format. "
                         "Format must be one of: {1}".format(file_format,
                                                             formats))


cdef class cDBG:

    SAVE_FORMATS = ['graphml', 'edgelist',
                    'gfa1', 'gfa2', 'fasta', 'gml']

    @staticmethod
    cdef cDBG _wrap(_cDBG * ptr):
        cdef cDBG cdbg = cDBG()
        cdbg._this = ptr
        return cdbg

    def unodes(self):
        cdef _cDBG.unode_iter_t it = deref(self._this).unodes_begin()
        cdef _UnitigNode * unode
        while(it != deref(self._this).unodes_end()):
            unode = deref(it).second.get()
            yield UnitigNode._wrap(unode)
            princ(it)


    def dnodes(self):
        cdef _cDBG.dnode_iter_t it = deref(self._this).dnodes_begin()
        cdef _DecisionNode * dnode
        while(it != deref(self._this).dnodes_end()):
            dnode = deref(it).second.get()
            yield DecisionNode._wrap(dnode)
            princ(it)

    def query_unode_hash(self, hash_t h):
        cdef _UnitigNode * _node = deref(self._this).query_unode_tag(h)
        if _node != NULL:
            return UnitigNode._wrap(_node)
        else:
            return None

    def query_unode_id(self, id_t node_id):
        cdef _UnitigNode * _node = \
            deref(self._this).query_unode_id(node_id)
        if _node != NULL:
            return UnitigNode._wrap(_node)
        else:
            return None

    def query_unode_end(self, hash_t h):
        cdef _UnitigNode * _node = \
            deref(self._this).query_unode_end(h)
        if _node != NULL:
            return UnitigNode._wrap(_node)
        else:
            return None
    
    def query_dnode(self, hash_t h):
        cdef _DecisionNode * _node = \
            deref(self._this).query_dnode(h)
        if _node != NULL:
            return DecisionNode._wrap(_node)
        else:
            return None

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

    def save(self, str filename, str file_format):
        if file_format is None:
            return
        else:
            deref(self._this).write(_bstring(filename),
                                    convert_format(file_format))


cdef class StreamingCompactor:

    def __cinit__(self, dBG__BitStorage__DefaultShifter graph):
        self._graph = graph._this.get()
        
        # TODO properly template this
        if type(self) is StreamingCompactor:
            self._sc_this = \
                make_unique[_StreamingCompactor[DefaultDBG]](self._graph)
            self.cdbg = cDBG._wrap(&(deref(self._sc_this).cdbg))
            self.Notifier = EventNotifier._wrap(<_EventNotifier*>self._sc_this.get())

    def wait_on_updates(self):
        deref(self._sc_this).wait_on_updates()

    def is_decision_kmer(self, str kmer):
        cdef string _kmer = _bstring(kmer)
        return deref(self._sc_this).is_decision_kmer(_kmer)

    def find_decision_kmers(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        cdef vector[uint32_t] positions
        cdef vector[hash_t] hashes
        cdef vector[NeighborBundle] neighbors

        deref(self._sc_this).find_decision_kmers(_sequence,
                                                 positions,
                                                 hashes,
                                                 neighbors)

        return positions, hashes

    def find_new_segments(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        cdef vector[vector[hash_t]] _segment_kmers
        cdef vector[string] _segment_seqs
        cdef deque[kmer_t] _decision_kmers
        cdef deque[NeighborBundle] _decision_neighbors

        deref(self._sc_this).find_new_segments(_sequence,
                                               _segment_kmers,
                                               _segment_seqs,
                                               _decision_kmers,
                                               _decision_neighbors)

        decision_kmers = []
        cdef int i = 0
        for i in range(_decision_kmers.size()):
            decision_kmers.append((_decision_kmers[i].hash,
                                   _decision_kmers[i].kmer))

        return _segment_seqs, decision_kmers
