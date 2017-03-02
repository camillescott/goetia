# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8

from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr

from lemon cimport SmartDigraph, INVALID


def _test_SmartDigraph_basic():

    cdef SmartDigraph * g = new SmartDigraph()
    cdef SmartDigraph.Node u = g.addNode()
    cdef SmartDigraph.Node v = g.addNode()
    cdef SmartDigraph.Arc e = g.addArc(u, v)

    cdef int cnt = 0
    cdef SmartDigraph.NodeIt n = SmartDigraph.NodeIt(deref(g))
    while n != SmartDigraph.Node(INVALID):
        preinc(cnt)
        preinc(n)
    assert cnt == 2
