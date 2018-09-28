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


cdef class CompactNodeView:

    @staticmethod
    cdef CompactNodeView _wrap(_CompactNode * ptr):
        cdef CompactNodeView node = CompactNodeView()
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
        return '<{ntype}: addr={addr} ID={id} L={L} sequence={seq}>'.format(ntype=self.__class__.__name__,
                                                                            id=self.node_id,
                                                                            L=len(self),
                                                                            seq=self.sequence,
                                                                            addr=str(<size_t>self._cn_this))
    def __repr__(self):
        self._check_ptr()
        return str(self)

    @property
    def node_id(self):
        self._check_ptr()
        return deref(self._cn_this).node_id

    def __eq__(self, other):
        return self.node_id == other.node_id


cdef class CompactNode:

    @property
    def sequence(self):
        return self.sequence

    @property
    def node_id(self):
        return self.node_id

    def __eq__(self, other):
        return self.node_id == other.node_id

    def __len__(self):
        return self.sequence.size()

    def __str__(self):
        return '<{ntype}: ID={id} L={L} sequence={seq}>'.format(ntype=self.__class__.__name__,
                                                                id=self.node_id,
                                                                L=len(self),
                                                                seq=self.sequence)


cdef class DecisionNodeView(CompactNodeView):

    @staticmethod
    cdef DecisionNodeView _wrap(_DecisionNode * ptr):
        cdef DecisionNodeView node = DecisionNodeView()
        node._dn_this = ptr
        node._cn_this = <_CompactNode*>ptr
        return node

    def clone(self):
        return DecisionNode._create(self)

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

    #def __str__(self):
    #    self._check_ptr()
    #    return deref(self._dn_this).repr()


cdef class DecisionNode(CompactNode):

    @staticmethod
    cdef DecisionNode _create(DecisionNodeView view):
        cdef DecisionNode node = DecisionNode()
        node.node_id = view.node_id
        node.sequence = deref(view._cn_this).sequence
        node.count = view.count
        return node

    @property
    def count(self):
        return self.count


cdef class UnitigNodeView(CompactNodeView):

    @staticmethod
    cdef UnitigNodeView _wrap(_UnitigNode * ptr):
        cdef UnitigNodeView node = UnitigNodeView()
        node._un_this = ptr
        node._cn_this = <_CompactNode*>ptr
        return node

    def clone(self):
        return UnitigNode._create(self)

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

    #def __str__(self):
    #    self._check_ptr()
    #    return deref(self._un_this).repr()


cdef class UnitigNode(CompactNode):

    @staticmethod
    cdef UnitigNode _create(UnitigNodeView view):
        cdef UnitigNode node = UnitigNode()
        node.node_id = view.node_id
        node.sequence = deref(view._cn_this).sequence
        node.left_end = view.left_end
        node.right_end = view.right_end
        node.meta = deref(view._un_this).meta()
        return node

    @property
    def meta(self):
        return _ustring(node_meta_repr(self.meta))

    @property
    def left_end(self):
        return self.left_end

    @property
    def right_end(self):
        return self.right_end


include "cdbg.tpl.pyx.pxi"
