# boink/minimizers.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libc.stdint cimport uint64_t, int64_t

from libcpp.memory cimport unique_ptr
from libcpp.utility cimport pair
from libcpp.vector cimport vector

cdef extern from "boink/minimizers.hh" namespace "boink":

    cdef cppclass _RollingMin "boink::RollingMin" [T]:
        _RollingMin(int64_t)

        void reset()
        const int64_t window_size()

        pair[T, int64_t] update(T)

    cdef cppclass _InteriorMinimizer "boink::InteriorMinimizer" [T](_RollingMin):
        vector[pair[T, int64_t]] get_minimizers()
        vector[T] get_minimizer_values()


cdef class InteriorMinimizer:
    cdef unique_ptr[_InteriorMinimizer[uint64_t]] _this
