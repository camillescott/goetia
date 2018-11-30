# boink/parsing.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libc.stdint cimport uint16_t, uint32_t, uint64_t, int64_t
from libcpp cimport bool
from libcpp.memory cimport unique_ptr
from libcpp.string cimport string


cdef extern from "boink/parsing/parsing.hh" namespace "boink::parsing" nogil:

    cdef cppclass _Sequence "boink::parsing::Read":
        string name
        string description
        string sequence
        string quality
        string cleaned_seq

        void reset()
        void set_cleaned_seq()

    cdef struct _SequenceBundle "boink::parsing::ReadBundle":
        bool has_left
        bool has_right
        _Sequence left
        _Sequence right

cdef extern from "boink/parsing/readers.hh" namespace "boink::parsing" nogil:
    cdef cppclass _SplitPairedReader "boink::parsing::SplitPairedReader" [ParserType=*]:
        _SplitPairedReader(const string&,
                           const string&,
                           uint32_t,
                           bool)

        _SplitPairedReader(const string&,
                           const string&)

        bool is_complete() except +ValueError
        _SequenceBundle next() except +ValueError


cdef class Sequence:
    cdef _Sequence _obj

    @staticmethod
    cdef Sequence _wrap(_Sequence cseq)


cdef class SplitPairedReader:

    cdef unique_ptr[_SplitPairedReader] _this
    cdef readonly int min_length
    cdef readonly bool force_name_match

