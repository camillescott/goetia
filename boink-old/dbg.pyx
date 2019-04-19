# boink/dbg.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

include "dbg.tpl.pyx.pxi"

from boink.hashing import UKHShifter

cdef class PdBG(dBG):

    def __cinit__(self, int K, int partition_K, *args, **kwargs):
        self.storage_type = '_PartitionedStorage'

        cdef vector[string] kmers = UKHShifter.get_kmers(K, partition_K)
        self._ukhs = make_shared[_UKHS](partition_K, kmers)
        if not self._this:
            self._this = make_shared[DefaultPdBG](K, partition_K, self._ukhs)
            self.allocated = True

    def hash(self, str kmer):
        return deref(self._this).hash(_bstring(kmer))

    def hashes(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef vector[hash_t] _hashes = deref(self._this).get_hashes(_sequence)
        cdef hash_t h
        for h in _hashes:
            yield h

    #def neighbors(self, str root):
    #    cdef bytes _root = _bstring(root)
    #    cdef pair[vector[kmer_t], vector[kmer_t]] result = deref(self._this).neighbor_kmers(_root)
    #    
    #    cdef list left = []
    #    cdef list right = []
    #    for i in range(result.first.size()):
    #        left.append((result.first[i].hash, result.first[i].kmer))
    #    for i in range(result.second.size()):
    #        right.append((result.second[i].hash, result.second[i].kmer))
    #
    #    return left, right

    def insert(self, str kmer):
        cdef string _kmer = _bstring(kmer)
        return deref(self._this).insert(_kmer)

    def insert_and_query(self, str kmer):
        cdef string _kmer = _bstring(kmer)
        return deref(self._this).insert_and_query(<string>_kmer)

    def query(self, str kmer):
        cdef string _kmer = _bstring(kmer)
        return deref(self._this).query(<string>_kmer)

    def insert_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        return deref(self._this).insert_sequence(<string>_sequence)       

    def insert_sequence_rolling(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        return deref(self._this).insert_sequence_rolling(<string>_sequence)   

    consume = insert_sequence

    def query_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef list counts = deref(self._this).query_sequence(_sequence)
        return counts

    def query_sequence_rolling(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef list counts = deref(self._this).query_sequence_rolling(_sequence)
        return counts

    @property
    def n_unique(self):
        return deref(self._this).n_unique()

    @property
    def n_occupied(self):
        return deref(self._this).n_occupied()

    @property
    def K(self):
        return deref(self._this).K()

    @property
    def partition_K(self):
        return deref(self._this).partition_K

    def save(self, file_name):
        deref(self._this).save(_bstring(file_name))

    def reset(self):
        deref(self._this).reset()

    def clone(self):
        cdef PdBG cloned = PdBG(self.K, self.partition_K)
        cloned._this = deref(self._this).clone()
        return cloned

    def get_partition_counts(self):
        return deref(self._this).get_partition_counts()
