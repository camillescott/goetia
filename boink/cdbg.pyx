# cdbg.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref
from libcpp.memory cimport make_shared

from boink.dbg cimport dBG_Base
from boink.utils cimport _bstring, _ustring
from khmer._oxli.sequence cimport Alphabets


cdef class CompactNode:

    @staticmethod
    cdef CompactNode _wrap(CompactNodePtr ptr):
        cdef CompactNode node = CompactNode()
        node._cn_this = ptr
        return node

    @property
    def sequence(self):
        return deref(self._cn_this).sequence

    def __len__(self):
        return deref(self._cn_this).sequence.length()

    def __str__(self):
        return 'CompactNode: L={0} sequence={1}'.format(len(self), self.sequence)

    def __repr__(self):
        return str(self)

    @property
    def node_id(self):
        return deref(self._cn_this).node_id


cdef class UnitigNode(CompactNode):

    @staticmethod
    cdef UnitigNode _wrap(UnitigNodePtr ptr):
        cdef UnitigNode node = UnitigNode()
        node._un_this = ptr
        return node

    @property
    def in_id(self):
        cdef uint64_t in_id = deref(self._un_this).in_id
        return None if in_id == NULL_ID else in_id

    @property
    def out_id(self):
        cdef uint64_t out_id = deref(self._un_this).out_id
        return None if out_id == NULL_ID else out_id

    def tags(self):
        cdef hash_t tag
        for tag in deref(self._un_this).tags:
            yield tag


cdef class DecisionNode(CompactNode):

    @staticmethod
    cdef DecisionNode _wrap(DecisionNodePtr ptr):
        cdef DecisionNode node = DecisionNode()
        node._dn_this = ptr
        return node

    @property
    def count(self):
        return deref(self._dn_this).count

    @property
    def out_degree(self):
        return deref(self._dn_this).out_degree()

    @property
    def in_degree(self):
        return deref(self._dn_this).in_degree()

    @property
    def degree(self):
        return deref(self._dn_this).degree()

    '''
    def out_edges(self):
        cdef string bases = Alphabets._get('DNA_SIMPLE')
        cdef char base
        cdef CpCompactEdge * edge
        for base in bases:
            edge = deref(self._cn_this).get_out_edge(base)
            if edge != NULL:
                yield <bytes>base, CompactEdge._wrap(edge)

    def in_edges(self):
        cdef string bases = Alphabets._get('DNA_SIMPLE')
        cdef char base
        cdef CpCompactEdge * edge
        for base in bases:
            edge = deref(self._cn_this).get_in_edge(base)
            if edge != NULL:
                yield <bytes>base, CompactEdge._wrap(edge)

    def __str__(self):
        return 'CompactNode: ID={0} count={1} in_degree={2}'\
               ' out_degree={3} sequence={4}'.format(self.kmer,
                                                     self.count,
                                                     self.in_degree,
                                                     self.out_degree,
                                                     self.sequence)
    '''


cdef class StreamingCompactor:

    def __cinit__(self, dBG_BitStorage_DefaultShifter graph):
        self._graph = graph._this
        
        # TODO properly template this
        if type(self) is StreamingCompactor:
            self._sc_this = \
                make_shared[_StreamingCompactor[DefaultDBG]](self._graph)

    def compactify(self, str seed):
        cdef bytes _seed = _bstring(seed)
        return deref(self._sc_this).compactify(_seed)


    def insert_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef vector[uint32_t] positions
        cdef vector[hash_t] hashes
        cdef vector[NeighborBundle] neighbors

        deref(self._sc_this).insert_sequence(_sequence,
                                             positions,
                                             hashes,
                                             neighbors)

        return positions, hashes

    '''
    def update(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        return deref(self._sc_this).update_compact_dbg(_sequence)

    def consume(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        return deref(self._sc_this).consume_sequence(_sequence)

    def consume_and_update(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        return deref(self._sc_this).consume_sequence_and_update(_sequence)

    def sequence_nodes(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        cdef vector[CpCompactNode*] nodes = \
                deref(self._sc_this).get_nodes(_sequence)
        cdef CpCompactNode* node
        for node in nodes:
            yield CompactNode._wrap(node)

    def report(self):
        return deref(self._sc_this).report()

    @property
    def n_nodes(self):
        return deref(self._sc_this).n_nodes()

    @property
    def n_edges(self):
        return deref(self._sc_this).n_edges()

    @property
    def n_updates(self):
        return deref(self._sc_this).n_updates()

    def write_gml(self, str filename):
        cdef string _filename = _bstring(filename)
        deref(self._sc_this).write_gml(_filename)

    def write_fasta(self, str filename):
        cdef string _filename = _bstring(filename)
        deref(self._sc_this).write_fasta(_filename)
    '''
