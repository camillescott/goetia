# boink/metrics.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libc.stdint   cimport uint64_t
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

cdef extern from "boink/metrics.hh" namespace "boink::metrics" nogil:
    cdef cppclass _ReservoirSample "boink::metrics::ReservoirSample" [T]:
        _ReservoirSample(size_t)

        void      sample(T)
        vector[T] get_result() const
        size_t    get_n_sampled() const
        size_t    get_sample_size() const
        void      clear()


cdef class ReservoirSample_Int:
    cdef unique_ptr[_ReservoirSample[uint64_t]] _this
