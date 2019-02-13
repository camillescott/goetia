# boink/hashing.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref

import os

from libc.stdint cimport uint64_t, UINT64_MAX
from libcpp.memory cimport unique_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector

from boink.utils cimport _bstring, _ustring
from boink.metadata import DATA_DIR


cdef void _test():
    cdef _DefaultShifter * shifter = new _DefaultShifter(4)
    print(deref(shifter).set_cursor(_bstring('AAAA')))
    print(deref(shifter).get_cursor())
    print(deref(shifter).shift_left(b'C'))
    print(deref(shifter).get_cursor())

    print(deref(shifter).hash(_bstring('AAAA')))
    print(deref(shifter).hash(_bstring('CAAA')))

def test():
    _test()


cdef class RollingHashShifter:

    cdef unique_ptr[_DefaultShifter] _this

    def __cinit__(self, uint16_t k):
        self._this.reset(new _DefaultShifter(k))

    def set_cursor(self, str kmer):
        deref(self._this).set_cursor(_bstring(kmer))

    def get_cursor(self):
        return deref(self._this).get_cursor()

    def hash(self, str kmer):
        return deref(self._this).hash(_bstring(kmer))

    def gather_left(self):
        cdef vector[shift_t] left = deref(self._this).gather_left()
        cdef shift_t item
        result = [(item.hash, _ustring(item.symbol)) for item in left]
        return result

    def gather_right(self):
        cdef vector[shift_t] right = deref(self._this).gather_right()
        cdef shift_t item
        result = [(item.hash, _ustring(item.symbol)) for item in right]
        return result

    @property
    def hashvalue(self):
        return deref(self._this).get()


cdef class UKHS:

    cdef unique_ptr[_UKHSShifter] _this

    def __cinit__(self, uint16_t W, uint16_t K):
        cdef vector[string] kmers = UKHS.get_kmers(W, K)
        self._this.reset(new _UKHSShifter(W, K, kmers))

    def set_cursor(self, str sequence):
        deref(self._this).set_cursor(_bstring(sequence))
        deref(self._this).reset_unikmer()

    @property
    def partition(self):
        return deref(self._this).get_partition()

    @property
    def unikmer(self):
        return deref(self._this).get_unikmer()

    @property
    def hashvalue(self):
        return deref(self._this).get()

    def query(self, uint64_t uhash):
        cdef uint64_t pos = UINT64_MAX
        cdef bool exists = deref(self._this).query(uhash, pos)
        if exists is True:
            return pos
        else:
            return None

    def hash(self, str kmer):
        return deref(self._this).hash(_bstring(kmer))

    @property
    def hashes(self):
        return deref(self._this).get_hashes()

    @staticmethod
    def get_kmers(int W, int K):
        valid_W = list(range(20, 210, 10))
        valid_K = list(range(7, 10))
        W = round(W, -1)
        if not W in valid_W:
            raise ValueError('Invalid UKHS window size.')
        if not K in valid_K:
            raise ValueError('Invalid UKHS K.')

        filename = os.path.join(DATA_DIR,
                                'res_{0}_{1}_4_0.txt'.format(K, W))
        kmers = []
        with open(filename) as fp:
            for line in fp:
                kmers.append(_bstring(line.strip()))

        return kmers
