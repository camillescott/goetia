# boink/dbg.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

cimport cython

from libc.stdint cimport uint8_t, uint16_t, uint64_t

from libcpp cimport bool
from libcpp.memory cimport shared_ptr, unique_ptr
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

cdef extern from "boink/storage/partitioned_storage.hh" nogil:
    cdef cppclass _PartitionedStorage "boink::storage::PartitionedStorage" [BaseStorage] (_Storage):
        pass

cdef extern from "boink/dbg.hh" namespace "boink" nogil:
    cdef cppclass _dBG "boink::dBG" [StorageType, HashShifter] (_KmerClient):
        _dBG(uint16_t)
        _dBG(uint16_t, vector[uint64_t])
        _dBG(uint16_t, uint64_t, uint16_t)
        _dBG(uint16_t, unique_ptr[StorageType]&)

        shared_ptr[_dBG[StorageType, HashShifter]] clone()

        ctypedef HashShifter shifter_type

        hash_t hash(string&) except +ValueError

        const bool insert(hash_t) 
        const bool insert(string&) except +ValueError
        const count_t insert_and_query(hash_t) 
        const count_t insert_and_query(string&) except +ValueError

        const count_t query(string&) except +ValueError
        const count_t query(hash_t)

        const string suffix(const string&)
        const string prefix(const string&)

        vector[kmer_t] left_neighbor_kmers(const string&)
        vector[kmer_t] right_neighbor_kmers(const string&)
        pair[vector[kmer_t], vector[kmer_t]] neighbor_kmers(const string&)

        vector[shift_t] left_neighbors(const string&)
        vector[shift_t] right_neighbors(const string&)
        pair[vector[shift_t], vector[shift_t]] neighbors(const string&)

        uint64_t insert_sequence(string&) except +ValueError
        uint64_t insert_sequence(string&,
                                 vector[hash_t]&,
                                 vector[count_t]&) except +ValueError

        vector[count_t] query_sequence(string&) except +ValueError
        void query_sequence(string&,
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

cdef extern from "boink/pdbg.hh" namespace "boink" nogil:
    cdef cppclass _PdBG "boink::PdBG" [BaseStorageType] (_KmerClient):
        _PdBG(uint16_t, uint16_t, shared_ptr[_UKHS], uint64_t, uint16_t)
        _PdBG(uint16_t, uint16_t, shared_ptr[_UKHS])

        const uint16_t partition_K

        shared_ptr[_PdBG[BaseStorageType]] clone() const

        const bool insert(const string&) except +ValueError
        const count_t insert_and_query(const string&) except +ValueError
        const count_t query(const string&) except +ValueError

        hash_t hash(string&) except +ValueError
        vector[hash_t] get_hashes(string&) except +ValueError

        const string suffix(const string&)
        const string prefix(const string&)

        vector[kmer_t] left_neighbor_kmers(const string&)
        vector[kmer_t] right_neighbor_kmers(const string&)
        pair[vector[kmer_t], vector[kmer_t]] neighbor_kmers(const string&)

        vector[shift_t] left_neighbors(const string&)
        vector[shift_t] right_neighbors(const string&)
        pair[vector[shift_t], vector[shift_t]] neighbors(const string&)

        uint64_t insert_sequence(string&) except +ValueError
        uint64_t insert_sequence(string&,
                                 vector[hash_t]&,
                                 vector[count_t]&) except +ValueError

        vector[count_t] query_sequence(string&) except +ValueError
        void query_sequence(string&,
                        vector[count_t]&,
                        vector[hash_t]&,
                        set[hash_t]&) except +ValueError

        uint64_t n_unique()
        uint64_t n_occupied()
        vector[size_t] get_partition_counts()

        void save(string)
        void load(string)
        void reset()

        shared_ptr[_KmerIterator[_DefaultUKHSShifter]] get_hash_iter(string&)

    ctypedef _PdBG[_SparseppSetStorage] DefaultPdBG
        


include "dbg.tpl.pxd.pxi"

cdef class PdBG(dBG):
    cdef shared_ptr[DefaultPdBG] _this
    cdef shared_ptr[_UKHS]       _ukhs

