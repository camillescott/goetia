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

from khmer._oxli.parsing cimport CpSequence as _Sequence


cdef extern from "boink/parsing.hh" namespace "boink" nogil:

    cdef struct _SequenceBundle "boink::ReadBundle":
        bool has_left
        bool has_right
        _Sequence left
        _Sequence right

    cdef cppclass _SplitPairedReader "boink::SplitPairedReader" [ParserType=*]:
        _SplitPairedReader(const string&,
                           const string&,
                           uint32_t,
                           bool)

        _SplitPairedReader(const string&,
                           const string&)

        bool is_complete() except +ValueError
        _SequenceBundle next() except +ValueError

cdef class SplitPairedReader:

    cdef unique_ptr[_SplitPairedReader] _this
    cdef readonly int min_length
    cdef readonly bool force_name_match

