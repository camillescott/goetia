# boink/minimizers.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libc.stdint cimport uint16_t, uint64_t, int64_t

from libcpp.memory cimport unique_ptr, shared_ptr
from libcpp.string cimport string
from libcpp.utility cimport pair
from libcpp.vector cimport vector

from boink.hashing cimport *
from boink.kmers cimport *


cdef extern from "boink/minimizers.hh" namespace "boink" nogil:

    cdef cppclass _RollingMin "boink::RollingMin" [T]:
        _RollingMin(int64_t)

        void reset()
        const int64_t window_size()
        const size_t size() const

        pair[T, int64_t] update(T)

    cdef cppclass _InteriorMinimizer "boink::InteriorMinimizer" [T](_RollingMin):
        vector[pair[T, int64_t]] get_minimizers()
        vector[T] get_minimizer_values()

    cdef cppclass _WKMinimizer "boink::WKMinimizer" [ShifterType](_InteriorMinimizer[hash_t], _KmerClient):
        _WKMinimizer(int64_t, uint16_t)

        vector[pair[hash_t, int64_t]] get_minimizers(const string&)
        vector[pair[string, int64_t]] get_minimizer_kmers(const string&)

cdef extern from "boink/ukhs_signature.hh" namespace "boink::signatures" nogil:

    cdef cppclass _UKHSSignature "boink::signatures::UKHSSignature" (_KmerClient):
        _UKHSSignature(uint16_t, uint16_t, uint64_t, shared_ptr[_UKHS])

        void insert(const string&) except +ValueError
        void insert_sequence(const string&) except +ValueError

        vector[vector[hash_t]] get_signature()
        uint64_t               get_n_accepted()
        uint64_t               get_n_rejected()
        uint64_t               get_n_kmers()
        vector[uint64_t]       get_bucket_n_inserts()
        

cdef class InteriorMinimizer:
    cdef unique_ptr[_InteriorMinimizer[hash_t]] _im_this
    cdef _InteriorMinimizer[hash_t] * _this


cdef class WKMinimizer(InteriorMinimizer):
    cdef unique_ptr[_WKMinimizer[_DefaultShifter]] _wk_this
