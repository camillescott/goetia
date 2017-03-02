# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8

from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr

from lemon cimport SmartDigraph, Node, NodeIt, Arc, INVALID


def _test_SmartDigraph_basic():

    cdef SmartDigraph * g = new SmartDigraph()
    cdef Node u = g.addNode()
    cdef Node v = g.addNode()
    cdef Arc e = g.addArc(u, v)

    cdef int cnt = 0
    cdef NodeIt n = NodeIt(deref(g))
    while n != INVALID:
        preinc(cnt)
        preinc(n)
    assert cnt == 2

