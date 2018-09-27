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


include "cdbg.tpl.pyx.pxi"
