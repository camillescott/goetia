# boink/dbg.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

cimport cython

from libc.stdint cimport uint8_t, uint16_t, uint64_t

from libcpp cimport bool
from libcpp.memory cimport shared_ptr
from libcpp.set cimport set
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

from boink.hashing cimport *

cdef extern from "boink/storage/storage.hh" namespace "boink::storage" nogil:
    # Need these for the Storage template parameter;
    # they don't need methods
    ctypedef uint8_t count_t
    ctypedef pair[count_t, count_t] full_count_t

    cdef cppclass _Storage "boink::storage::Storage":
        pass

cdef extern from "boink/storage/bitstorage.hh" nogil:
    cdef cppclass _BitStorage "boink::storage::BitStorage" (_Storage):
        pass

cdef extern from "boink/storage/nibblestorage.hh" nogil:
    cdef cppclass _NibbleStorage "boink::storage::NibbleStorage" (_Storage):
        pass

cdef extern from "boink/storage/qfstorage.hh" nogil:
    cdef cppclass _QFStorage "boink::storage::QFStorage" (_Storage):
        pass

cdef extern from "boink/storage/bytestorage.hh" nogil:
    cdef cppclass _ByteStorage "boink::storage::ByteStorage" (_Storage):
        pass

cdef extern from "boink/storage/sparseppstorage.hh" nogil:
    cdef cppclass _SparseppSetStorage "boink::storage::SparseppSetStorage" (_Storage):
        pass


cdef extern from "boink/dbg.hh" namespace "boink" nogil:
    cdef cppclass _dBG "boink::dBG" [StorageType, HashShifter] (_KmerClient):
        _dBG(uint16_t, vector[uint64_t])
        _dBG(uint16_t, uint64_t, uint16_t)

        shared_ptr[_dBG[StorageType, HashShifter]] clone()

        ctypedef HashShifter shifter_type

        hash_t hash(string&) except +ValueError
        bool add(hash_t)
        bool add(string&) except +ValueError

        count_t get(string&) except +ValueError
        count_t get(hash_t)

        const string suffix(const string&)
        const string prefix(const string&)

        vector[kmer_t] left_neighbor_kmers(const string&)
        vector[kmer_t] right_neighbor_kmers(const string&)
        pair[vector[kmer_t], vector[kmer_t]] neighbor_kmers(const string&)

        vector[shift_t] left_neighbors(const string&)
        vector[shift_t] right_neighbors(const string&)
        pair[vector[shift_t], vector[shift_t]] neighbors(const string&)

        uint64_t add_sequence(string&,
                              vector[hash_t]&,
                              vector[bool]&) except +ValueError
        uint64_t add_sequence(string&) except +ValueError
        vector[count_t] get_counts(string&) except +ValueError
        void get_counts(string&,
                        vector[count_t]&,
                        vector[hash_t]&,
                        set[hash_t]&) except +ValueError

        uint64_t n_unique()
        uint64_t n_occupied()

        uint8_t ** get_raw()

        void save(string)
        void load(string)
        void reset()

        shared_ptr[_KmerIterator[HashShifter]] get_hash_iter(string&)

    ctypedef _dBG[_BitStorage, _DefaultShifter] DefaultDBG


include "dbg.tpl.pxd.pxi"
