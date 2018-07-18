# boink/assembly.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libc.stdint cimport uint8_t

from libcpp cimport bool
from libcpp.deque cimport deque
from libcpp.string cimport string
from libcpp.vector import vector

from boink.hashing cimport *


cdef extern from "boink/assembly.hh" namespace "boink" nogil:
    ctypedef deque[char] Path
    ctypedef vector[string] StringVector
    ctypedef vector[kmer_t] KmerVector

    cdef cppclass _AssemblerMixin "boink::AssemblerMixin" [GraphType]:

        _AssemblerMixin(GraphType *)

        void clear_seen()

        uint8_t degree_left()
        uint8_t degree_right()
        uint8_t degree()

        bool get_left(shift_t&)
        bool get_right(shift_t&)
        bool reduce_nodes(vector[shift_t]&, shift_t&)
        uint8_t count_nodes(vector[shift_t]&)
        vector[shift_t] filter_nodes(vector[shift_t]&)

        void assemble_left(const string&, Path&) except +ValueError
        void assemble_left(Path&)

        void assemble_right(const string&, Path&) except +ValueError
        void assemble_right(Path&)

        void assemble(const string&, Path&) except +ValueError

        string to_string(Path&)

        # HashShifter methods
        hash_t set_cursor(string&) except +ValueError
        string get_cursor()
        void get_cursor(deque[char]&)

        bool is_valid(const char)
        bool is_valid(const string&)

        hash_t get()
        hash_t hash(string&) except +ValueError

        vector[shift_t] gather_left()
        vector[shift_t] gather_right()

        hash_t shift_left(const char) except +ValueError
        hash_t shift_right(const char) except +ValueError

    _AssemblerMixin[GraphType] make_assembler[GraphType](GraphType *)


include "assembly.tpl.pxd.pxi"
