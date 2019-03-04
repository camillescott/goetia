# boink/minimizers.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref

from libcpp.memory cimport make_unique, make_shared

from boink.utils cimport _bstring
from boink.hashing import UKHShifter

cdef class InteriorMinimizer:

    def __cinit__(self, int64_t window_size, *args, **kwargs):
        if type(self) is InteriorMinimizer:
            self._im_this = make_unique[_InteriorMinimizer[hash_t]](window_size)
            self._this = self._im_this.get()

    @property
    def window_size(self):
        return deref(self._this).window_size()

    def reset(self):
        deref(self._this).reset()

    def update(self, hash_t value):
        '''Update the minimizer with a new value. Returns 
        tuple of (current_min, index).'''

        cdef pair[hash_t, int64_t] result = deref(self._this).update(value)
        return result.first, result.second

    def __call__(self, object values):
        for value in values:
            self.update(value)
        return self.get_minimizer_values()

    def get_minimizers(self):
        cdef vector[pair[hash_t, int64_t]] minimizers = deref(self._this).get_minimizers()
        return minimizers

    def get_minimizer_values(self):
        cdef vector[hash_t] values = deref(self._this).get_minimizer_values()
        return values


cdef class WKMinimizer(InteriorMinimizer):
    
    def __cinit__(self, int64_t window_size, uint16_t ksize, *args, **kwargs):
        if type(self) is WKMinimizer:
            self._wk_this = make_unique[_WKMinimizer[_DefaultShifter]](window_size, ksize)
            self._this = <_InteriorMinimizer[hash_t]*>self._wk_this.get()

    def get_minimizers(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        return deref(self._wk_this).get_minimizers(_sequence)

    def get_minimizer_kmers(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        return deref(self._wk_this).get_minimizer_kmers(_sequence)


cdef class UKHSSignature:

    cdef shared_ptr[_UKHSSignature] _this
    cdef shared_ptr[_UKHS]          _ukhs
    cdef readonly int               K
    cdef readonly int               bucket_K
    cdef readonly int               bucket_size

    def __cinit__(self, int K, int bucket_K, int bucket_size):
        cdef vector[string] kmers = UKHShifter.get_kmers(K, bucket_K)
        self._ukhs = make_shared[_UKHS](bucket_K, kmers)
        if not self._this:
            self._this       = make_shared[_UKHSSignature](K, bucket_K, bucket_size, self._ukhs)
            self.K           = K
            self.bucket_K    = bucket_K
            self.bucket_size = bucket_size
    
    def insert_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        return deref(self._this).insert_sequence(_sequence)       


    @property
    def signature(self):
        sig = deref(self._this).get_signature()
        return sig

