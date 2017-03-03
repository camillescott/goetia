# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8, embedsignature=True

from dbg cimport ExactDBG
from dbg import ExactDBG

from khmer._oxli.parsing cimport SanitizedFastxParser as SFParser, Sequence


cdef class FileConsumer:

    def __cinit__(self, ExactDBG graph):
        self.graph = graph

    cpdef int consume(self, filename) except -1:
        raise NotImplementedError()


cdef class FastxConsumer(FileConsumer):

    def __cinit__(self, ExactDBG graph):
        self.graph = graph

    cpdef int consume(self, filename) except -1:
        cdef SFParser parser = SFParser(filename, 
                                        alphabet=self.graph.alphabet,
                                        convert_n=False)
        cdef Sequence seq
        cdef int n_consumed = 0
        while not parser.is_complete():
            seq = parser._next()
            if seq is not None:
                self.graph._add_sequence(seq._obj.sequence)
                n_consumed += 1
        
        return n_consumed
