# boink/metrics.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref
from libcpp.memory cimport make_unique


cdef class ReservoirSample_Int:

    def __cinit__(self, size_t sample_size):
        if type(self) is ReservoirSample_Int:
            self._this = make_unique[_ReservoirSample[uint64_t]](sample_size)

    @property
    def n_sampled(self):
        return deref(self._this).get_n_sampled()

    @property
    def sample_size(self):
        return deref(self._this).get_sample_size()

    def clear(self):
        deref(self._this).clear()

    def sample(self, uint64_t item):
        deref(self._this).sample(item)

    def get_result(self):
        return deref(self._this).get_result()
