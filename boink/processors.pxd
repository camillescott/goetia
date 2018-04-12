# boink/processors.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libc.stdint cimport uint16_t, uint32_t, uint64_t, int64_t
from libcpp.string cimport string

from boink.dbg cimport *
from boink.cdbg cimport *
from boink.utils cimport _bstring

cdef extern from "boink/consumer.hh" namespace "boink":

    cdef cppclass _FileProcessor "boink::FileProcessor" [Derived]:
        void process(const string&,
                     uint64_t&,
                     uint64_t&,
                     uint32_t) except +ValueError
        void process(const string&,
                     uint64_t&,
                     uint64_t&) except +ValueError

    cdef cppclass _FileConsumer "boink::FileConsumer" [GraphType] (_FileProcessor[_FileConsumer[GraphType]]):
        _FileConsumer(GraphType *)

    cdef cppclass _DecisionNodeProcessor "boink::DecisionNodeProcessor"[GraphType] (_FileProcessor[_DecisionNodeProcessor[GraphType]]):
        _DecisionNodeProcessor(_StreamingCompactor*,
                               string &)

    cdef cppclass _MinimizerProcessor "boink::MinimizerProcessor" [ShifterType] (_FileProcessor[_MinimizerProcessor[ShifterType]]):
        _MinimizerProcessor(int64_t, uint16_t, const string&)


cdef class FileProcessor:
    pass

cdef class FileConsumer(FileProcessor):
    cdef unique_ptr[_FileConsumer[DefaultDBG]] _fc_this


cdef class DecisionNodeProcessor(FileProcessor):
    cdef readonly str output_filename
    cdef unique_ptr[_DecisionNodeProcessor[DefaultDBG]] _dnp_this


cdef class MinimizerProcessor(FileProcessor):
    cdef unique_ptr[_MinimizerProcessor[DefaultShifter]] _mp_this
