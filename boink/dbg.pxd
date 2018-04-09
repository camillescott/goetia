# boink/dbg.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

cimport cython

from libc.stdint cimport uint8_t, uint16_t, uint64_t

from libcpp cimport bool
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

from boink.hashing cimport *

cdef extern from "oxli/storage.hh":
    # Need these for the Storage template parameter;
    # they don't need methods
    cdef cppclass Storage "oxli::Storage":
        pass
    cdef cppclass BitStorage "oxli::BitStorage" (Storage):
        pass
    cdef cppclass NibbleStorage "oxli::NibbleStorage" (Storage):
        pass
    cdef cppclass QFStorage "oxli::QFStorage" (Storage):
        pass
    cdef cppclass ByteStorage "oxli::ByteStorage" (Storage):
        pass


cdef extern from "boink/dbg.hh" namespace "boink":
    ctypedef pair[bool, bool] bit_pair_t
    ctypedef vector[bit_pair_t] bit_pair_vector_t
    ctypedef uint8_t count_t
    ctypedef pair[count_t, count_t] full_count_t

    cdef cppclass _dBG "boink::dBG" [StorageType, HashShifter] (_KmerClient):
        _dBG(uint16_t, vector[uint64_t])
        _dBG(uint16_t, uint64_t, uint16_t)

        unique_ptr[_dBG[StorageType, HashShifter]] clone()

        ctypedef HashShifter shifter_type

        hash_t hash(string&) except +ValueError
        bool add(hash_t)
        bool add(string&) except +ValueError

        count_t get(string&) except +ValueError
        count_t get(hash_t)

        uint64_t add_sequence(string&, vector[bool]&) except +ValueError
        uint64_t add_sequence(string&) except +ValueError
        vector[count_t] get_counts(string&) except +ValueError

        uint64_t n_unique()
        uint64_t n_occupied()

        uint8_t ** get_raw()

        void save(string)
        void load(string)
        void reset()

        unique_ptr[_KmerIterator[HashShifter]] get_hash_iter(string&)

    ctypedef _dBG[BitStorage, DefaultShifter] DefaultDBG


include "dbg.pxd.pxi"
