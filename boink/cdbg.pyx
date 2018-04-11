# cdbg.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref
from libcpp.memory cimport make_unique

from boink.utils cimport _bstring, _ustring
from khmer._oxli.sequence cimport Alphabets


cdef class CompactNode:

    @staticmethod
    cdef CompactNode _wrap(_CompactNode * ptr):
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


cdef class DecisionNode(CompactNode):

    @staticmethod
    cdef DecisionNode _wrap(_DecisionNode * ptr):
        cdef DecisionNode node = DecisionNode()
        node._dn_this = ptr
        node._cn_this = <_CompactNode*>ptr
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

    def in_edges(self):
        cdef id_t node_id
        for node_id in deref(self._dn_this).in_edges:
            yield node_id


    def out_edges(self):
        cdef id_t node_id
        for node_id in deref(self._dn_this).out_edges:
            yield node_id

    '''


    def __str__(self):
        return 'CompactNode: ID={0} count={1}'
               ' in_degree={2}'
               ' out_degree={3} sequence={4}'.format(self.kmer,
                                                     self.count,
                                                     self.in_degree,
                                                     self.out_degree,
                                                     self.sequence)
    '''

cdef class UnitigNode(CompactNode):

    @staticmethod
    cdef UnitigNode _wrap(_UnitigNode * ptr):
        cdef UnitigNode node = UnitigNode()
        node._un_this = ptr
        node._cn_this = <_CompactNode*>ptr
        return node

    @property
    def left_hash(self):
        return deref(self._un_this).left_hash

    @property
    def right_hash(self):
        return deref(self._un_this).right_hash

    @property
    def left_dnode(self):
        cdef _DecisionNode * ld = deref(self._un_this).left_dnode
        return None if ld == NULL else DecisionNode._wrap(ld)

    @property
    def right_dnode(self):
        cdef _DecisionNode * rd = deref(self._un_this).right_dnode
        return None if rd == NULL else DecisionNode._wrap(rd)

    def tags(self):
        cdef hash_t tag
        for tag in deref(self._un_this).tags:
            yield tag


cdef class StreamingCompactor:

    def __cinit__(self, dBG_BitStorage_DefaultShifter graph):
        self._graph = graph._this.get()
        
        # TODO properly template this
        if type(self) is StreamingCompactor:
            self._sc_this = \
                make_unique[_StreamingCompactor[DefaultDBG]](self._graph)
            self._cdbg = &(deref(self._sc_this).cdbg)

    #def compactify(self, str seed):
    #    cdef string _seed = _bstring(seed)
    #    return deref(self._sc_this).compactify(_seed)

    def is_decision_node(self, str kmer):
        cdef string _kmer = _bstring(kmer)
        return deref(self._sc_this).is_decision_node(_kmer)

    def find_decision_nodes(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        cdef vector[uint32_t] positions
        cdef vector[hash_t] hashes
        cdef vector[NeighborBundle] neighbors

        deref(self._sc_this).find_decision_nodes(_sequence,
                                                 positions,
                                                 hashes,
                                                 neighbors)

        return positions, hashes

    def get_cdbg_dnodes(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        cdef vector[_DecisionNode*] dnodes = \
            deref(self._cdbg).get_decision_nodes[DefaultShifter](_sequence)
        cdef _DecisionNode * dnode
        for dnode in dnodes:
            if dnode != NULL:
                yield DecisionNode._wrap(dnode)

    def get_cdbg_unode_from_id(self, id_t node_id):
        cdef _UnitigNode * unode = deref(self._cdbg).get_unitig_node_from_id(node_id)
        return None if unode == NULL else UnitigNode._wrap(unode)
        

    def insert_sequence(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        cdef vector[uint32_t] positions
        cdef vector[hash_t] hashes
        cdef vector[NeighborBundle] neighbors

        deref(self._sc_this).insert_sequence(_sequence,
                                             positions,
                                             hashes,
                                             neighbors)

        return positions, hashes

    def update(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        return deref(self._sc_this).update(_sequence)

    '''
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
