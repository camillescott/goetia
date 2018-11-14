# boink/utils.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libcpp cimport bool
from libcpp.utility cimport pair


cdef extern from "<utility>" namespace "std" nogil:
    pair[U,T] make_pair[U,T](...) except +

cdef bytes _bstring(s)
cdef unicode _ustring(s)

cpdef bool is_str(object s)
cpdef bool is_num(object n)


