# boink/utils.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libcpp cimport bool
from libcpp.utility cimport pair
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint32_t, uint64_t

cdef extern from "<utility>" namespace "std" nogil:
    pair[U,T] make_pair[U,T](...) except +

cdef bytes _bstring(s)
cdef unicode _ustring(s)

cpdef bool is_str(object s)
cpdef bool is_num(object n)

cpdef bool is_prime(uint64_t n)
cdef vector[uint64_t] get_n_primes_near_x(uint32_t, uint64_t)
