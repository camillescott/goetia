# boink/kmers.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libc.stdint cimport uint16_t

cdef extern from "boink/kmers/kmerclient.hh" namespace "boink::kmers" nogil:

    cdef cppclass _KmerClient "boink::kmers::KmerClient":
        _KmerClient(uint16_t)
        uint16_t K()
