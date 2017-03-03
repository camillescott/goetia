# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8, embedsignature=True

from dbg cimport ExactDBG


cdef class FileConsumer:

    cdef ExactDBG graph

    cpdef int consume(self, filename) except -1

cdef class FastxConsumer(FileConsumer):
    pass
