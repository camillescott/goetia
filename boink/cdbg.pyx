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
        return deref(self._dn_this).count

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

    def left_junctions(self):
        self._check_ptr()
        cdef junction_t l_junc
        for l_junc in deref(self._dn_this).left_juncs:
            yield l_junc.first, l_junc.second


    def right_junctions(self):
        self._check_ptr()
        cdef junction_t r_junc
        for r_junc in deref(self._dn_this).right_juncs:
            yield r_junc.first, r_junc.second

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
    def left_junc(self):
        self._check_ptr()
        return (deref(self._un_this).left_junc.first,
                deref(self._un_this).left_junc.second)

    @property
    def right_junc(self):
        self._check_ptr()
        return (deref(self._un_this).right_junc.first,
                deref(self._un_this).right_junc.second)

    @property
    def meta(self):
        return _ustring(node_meta_repr(deref(self._un_this).meta))

    def tags(self):
        self._check_ptr()
        cdef hash_t tag
        for tag in deref(self._un_this).tags:
            yield tag

    def __str__(self):
        self._check_ptr()
        return deref(self._un_this).repr()


cdef class cDBG:

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

    def get_unode_by_hash(self, hash_t h):
        cdef _UnitigNode * _node = deref(self._this).get_unode(h)
        if _node != NULL:
            return UnitigNode._wrap(_node)
        else:
            return None

    def get_unode_by_id(self, id_t node_id):
        cdef _UnitigNode * _node = \
            deref(self._this).get_unode_from_id(node_id)
        if _node != NULL:
            return UnitigNode._wrap(_node)
        else:
            return None

    def get_unode_by_junc(self, tuple junc):
        cdef junction_t j = make_pair[hash_t, hash_t](<hash_t>junc[0],
                                                      <hash_t>junc[1])
        cdef _UnitigNode * _node = deref(self._this).get_unode(j)
        if _node != NULL:
            return UnitigNode._wrap(_node)
        else:
            return None

    def get_dnode(self, hash_t h):
        cdef _DecisionNode * _node = \
            deref(self._this).get_dnode(h)
        if _node != NULL:
            return DecisionNode._wrap(_node)
        else:
            return None

    def get_dnodes(self, str sequence):
        pass

    def get_left_dnode(self, UnitigNode unode):
        cdef _DecisionNode * dnode = deref(self._this).get_left_dnode(
                                                unode._un_this
                                                )
        return None if dnode == NULL else DecisionNode._wrap(dnode)

    def get_right_dnode(self, UnitigNode unode):
        cdef _DecisionNode * dnode = deref(self._this).get_right_dnode(
                                                unode._un_this
                                                )
        return None if dnode == NULL else DecisionNode._wrap(dnode)

    def neighbors(self, CompactNode node):
        for neighbor in node.left_neighbors(node):
            yield neighbor
        for neighbor in self.right_neighbor(node):
            yield neighbor

    def left_neighbors(self, CompactNode node):
        if isinstance(node, UnitigNode):
            neighbor = self.get_left_dnode(node)
            if neighbor is not None:
                yield neighbor
        if isinstance(node, DecisionNode):
            for junc in node.left_junctions():
                yield self.get_unode_by_junc(junc)

    def right_neighbors(self, CompactNode node):
        if isinstance(node, UnitigNode):
            neighbor = self.get_right_dnode(node)
            if neighbor is not None:
                yield neighbor
        if isinstance(node, DecisionNode):
            for junc in node.right_junctions():
                yield self.get_unode_by_junc(junc)

    @property
    def n_updates(self):
        return deref(self._this).n_updates()

    @property
    def n_unodes(self):
        return deref(self._this).n_unitig_nodes()

    @property
    def n_dnodes(self):
        return deref(self._this).n_decision_nodes()


cdef class StreamingCompactor:

    def __cinit__(self, dBG__BitStorage__DefaultShifter graph):
        self._graph = graph._this.get()
        
        # TODO properly template this
        if type(self) is StreamingCompactor:
            self._sc_this = \
                make_unique[_StreamingCompactor[DefaultDBG]](self._graph)
            self.cdbg = cDBG._wrap(&(deref(self._sc_this).cdbg))

    #def compactify(self, str seed):
    #    cdef string _seed = _bstring(seed)
    #    return deref(self._sc_this).compactify(_seed)

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

    def get_cdbg_dnodes(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        cdef vector[_DecisionNode*] dnodes = \
            deref(self.cdbg._this).get_dnodes[_DefaultShifter](_sequence)
        cdef _DecisionNode * dnode
        for dnode in dnodes:
            if dnode != NULL:
                yield DecisionNode._wrap(dnode)

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

    def update_sequence(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        return deref(self._sc_this).update_sequence(_sequence)

    def update_cdbg(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        return deref(self._sc_this).update_cdbg(_sequence)

    '''
    def write_gml(self, str filename):
        cdef string _filename = _bstring(filename)
        deref(self._sc_this).write_gml(_filename)

    def write_fasta(self, str filename):
        cdef string _filename = _bstring(filename)
        deref(self._sc_this).write_fasta(_filename)
    '''
