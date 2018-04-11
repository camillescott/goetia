# boink/minimizers.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref

from libcpp.memory cimport make_unique

cdef class InteriorMinimizer:

    def __cinit__(self, int64_t window_size):
        self._this = make_unique[_InteriorMinimizer[uint64_t]](window_size)

    @property
    def window_size(self):
        return deref(self._this).window_size()

    def reset(self):
        deref(self._this).reset()

    def update(self, uint64_t value):
        '''Update the minimizer with a new value. Returns 
        tuple of (current_min, index).'''

        cdef pair[uint64_t, int64_t] result = deref(self._this).update(value)
        return result.first, result.second

    def __call__(self, object values):
        for value in values:
            self.update(value)
        return self.get_minimizer_values()

    def get_minimizers(self):
        cdef vector[pair[uint64_t, int64_t]] minimizers = deref(self._this).get_minimizers()
        return minimizers

    def get_minimizer_values(self):
        cdef vector[uint64_t] values = deref(self._this).get_minimizer_values()
        return values
