# boink/dbg.pyx.pxi
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.
#
# WARNING: this file is automatically generated; do not modify it!
# The source template is: dbg.tpl.pyx

from cython.operator cimport dereference as deref

from libc.stdint cimport uint64_t
from libcpp.string cimport string
from libcpp.vector cimport vector

from khmer._oxli.utils import get_n_primes_near_x, is_str, is_num
from boink.utils cimport _bstring, _ustring


cdef class dBG_Base:
    
    def __cinit__(self, *args, **kwargs):
        self.allocated = False

cdef class dBG_BitStorage_DefaultShifter(dBG_Base):

    def __cinit__(self, int K, uint64_t starting_size, int n_tables,
                  *args, **kwargs):
        cdef vector[uint64_t] primes
        #if type(self) is dBG_BitStorage_DefaultShifter:
        if not self._this:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self._this = make_shared[_dBG[BitStorage,DefaultShifter]](K, primes)
            self._assembler = make_shared[_AssemblerMixin[_dBG[BitStorage,DefaultShifter]]](self._this)
            self.allocated = True

        self.storage_type = "BitStorage"
        self.shifter_type = "DefaultShifter"
        self.suffix = "BitStorage_DefaultShifter"

    cdef hash_t _handle_kmer(self, object kmer) except 0:
        cdef hash_t handled
        if is_num(kmer):
            handled = <hash_t> kmer
        else:
            handled = deref(self._this).hash(_bstring(kmer))
        return handled

    def add(self, object kmer):
        return deref(self._this).add(self._handle_kmer(kmer))

    def get(self, object kmer):
        return deref(self._this).get(self._handle_kmer(kmer))

    def hash(self, str kmer):
        return deref(self._this).hash(_bstring(kmer))

    def hashes(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef shared_ptr[_KmerIterator[DefaultShifter]] kmer_iter = \
                deref(self._this).get_hash_iter(_sequence)

        cdef hash_t h
        while(not deref(kmer_iter).done()):
            h = deref(kmer_iter).next()
            yield h

    def add_sequence_and_report(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef vector[bool] report
        deref(self._this).add_sequence(_sequence, report)
        return report

    def add_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        return deref(self._this).add_sequence(_sequence)

    # compatibility with oxli API
    consume = add_sequence

    def get_counts(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef list counts = deref(self._this).get_counts(_sequence)
        return counts

    def left_degree(self, str kmer):
        deref(self._assembler).set_cursor(_bstring(kmer))
        return deref(self._assembler).degree_left()

    def right_degree(self, str kmer):
        deref(self._assembler).set_cursor(_bstring(kmer))
        return deref(self._assembler).degree_right()

    def degree(self, str kmer):
        deref(self._assembler).set_cursor(_bstring(kmer))
        return deref(self._assembler).degree_left() + \
               deref(self._assembler).degree_right()

    @property
    def n_unique(self):
        return deref(self._this).n_unique()

    @property
    def n_occupied(self):
        return deref(self._this).n_occupied()

    @property
    def K(self):
        return deref(self._this).K()

    def save(self, file_name):
        deref(self._this).save(_bstring(file_name))

    @classmethod
    def load(cls, file_name):
        cdef dBG_BitStorage_DefaultShifter obj = cls(1, 1, 1)
        deref(obj._this).load(_bstring(file_name))
        return obj

cdef class dBG_NibbleStorage_DefaultShifter(dBG_Base):

    def __cinit__(self, int K, uint64_t starting_size, int n_tables,
                  *args, **kwargs):
        cdef vector[uint64_t] primes
        #if type(self) is dBG_NibbleStorage_DefaultShifter:
        if not self._this:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self._this = make_shared[_dBG[NibbleStorage,DefaultShifter]](K, primes)
            self._assembler = make_shared[_AssemblerMixin[_dBG[NibbleStorage,DefaultShifter]]](self._this)
            self.allocated = True

        self.storage_type = "NibbleStorage"
        self.shifter_type = "DefaultShifter"
        self.suffix = "NibbleStorage_DefaultShifter"

    cdef hash_t _handle_kmer(self, object kmer) except 0:
        cdef hash_t handled
        if is_num(kmer):
            handled = <hash_t> kmer
        else:
            handled = deref(self._this).hash(_bstring(kmer))
        return handled

    def add(self, object kmer):
        return deref(self._this).add(self._handle_kmer(kmer))

    def get(self, object kmer):
        return deref(self._this).get(self._handle_kmer(kmer))

    def hash(self, str kmer):
        return deref(self._this).hash(_bstring(kmer))

    def hashes(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef shared_ptr[_KmerIterator[DefaultShifter]] kmer_iter = \
                deref(self._this).get_hash_iter(_sequence)

        cdef hash_t h
        while(not deref(kmer_iter).done()):
            h = deref(kmer_iter).next()
            yield h

    def add_sequence_and_report(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef vector[bool] report
        deref(self._this).add_sequence(_sequence, report)
        return report

    def add_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        return deref(self._this).add_sequence(_sequence)

    # compatibility with oxli API
    consume = add_sequence

    def get_counts(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef list counts = deref(self._this).get_counts(_sequence)
        return counts

    def left_degree(self, str kmer):
        deref(self._assembler).set_cursor(_bstring(kmer))
        return deref(self._assembler).degree_left()

    def right_degree(self, str kmer):
        deref(self._assembler).set_cursor(_bstring(kmer))
        return deref(self._assembler).degree_right()

    def degree(self, str kmer):
        deref(self._assembler).set_cursor(_bstring(kmer))
        return deref(self._assembler).degree_left() + \
               deref(self._assembler).degree_right()

    @property
    def n_unique(self):
        return deref(self._this).n_unique()

    @property
    def n_occupied(self):
        return deref(self._this).n_occupied()

    @property
    def K(self):
        return deref(self._this).K()

    def save(self, file_name):
        deref(self._this).save(_bstring(file_name))

    @classmethod
    def load(cls, file_name):
        cdef dBG_NibbleStorage_DefaultShifter obj = cls(1, 1, 1)
        deref(obj._this).load(_bstring(file_name))
        return obj

cdef class dBG_ByteStorage_DefaultShifter(dBG_Base):

    def __cinit__(self, int K, uint64_t starting_size, int n_tables,
                  *args, **kwargs):
        cdef vector[uint64_t] primes
        #if type(self) is dBG_ByteStorage_DefaultShifter:
        if not self._this:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self._this = make_shared[_dBG[ByteStorage,DefaultShifter]](K, primes)
            self._assembler = make_shared[_AssemblerMixin[_dBG[ByteStorage,DefaultShifter]]](self._this)
            self.allocated = True

        self.storage_type = "ByteStorage"
        self.shifter_type = "DefaultShifter"
        self.suffix = "ByteStorage_DefaultShifter"

    cdef hash_t _handle_kmer(self, object kmer) except 0:
        cdef hash_t handled
        if is_num(kmer):
            handled = <hash_t> kmer
        else:
            handled = deref(self._this).hash(_bstring(kmer))
        return handled

    def add(self, object kmer):
        return deref(self._this).add(self._handle_kmer(kmer))

    def get(self, object kmer):
        return deref(self._this).get(self._handle_kmer(kmer))

    def hash(self, str kmer):
        return deref(self._this).hash(_bstring(kmer))

    def hashes(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef shared_ptr[_KmerIterator[DefaultShifter]] kmer_iter = \
                deref(self._this).get_hash_iter(_sequence)

        cdef hash_t h
        while(not deref(kmer_iter).done()):
            h = deref(kmer_iter).next()
            yield h

    def add_sequence_and_report(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef vector[bool] report
        deref(self._this).add_sequence(_sequence, report)
        return report

    def add_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        return deref(self._this).add_sequence(_sequence)

    # compatibility with oxli API
    consume = add_sequence

    def get_counts(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef list counts = deref(self._this).get_counts(_sequence)
        return counts

    def left_degree(self, str kmer):
        deref(self._assembler).set_cursor(_bstring(kmer))
        return deref(self._assembler).degree_left()

    def right_degree(self, str kmer):
        deref(self._assembler).set_cursor(_bstring(kmer))
        return deref(self._assembler).degree_right()

    def degree(self, str kmer):
        deref(self._assembler).set_cursor(_bstring(kmer))
        return deref(self._assembler).degree_left() + \
               deref(self._assembler).degree_right()

    @property
    def n_unique(self):
        return deref(self._this).n_unique()

    @property
    def n_occupied(self):
        return deref(self._this).n_occupied()

    @property
    def K(self):
        return deref(self._this).K()

    def save(self, file_name):
        deref(self._this).save(_bstring(file_name))

    @classmethod
    def load(cls, file_name):
        cdef dBG_ByteStorage_DefaultShifter obj = cls(1, 1, 1)
        deref(obj._this).load(_bstring(file_name))
        return obj


cdef object _make_dbg(int K, uint64_t starting_size, int n_tables,
                      str storage='BitStorage',
                      str shifter='DefaultShifter'):

    if storage == "BitStorage" and shifter == "DefaultShifter":
        return dBG_BitStorage_DefaultShifter(K, starting_size, n_tables)

    if storage == "NibbleStorage" and shifter == "DefaultShifter":
        return dBG_NibbleStorage_DefaultShifter(K, starting_size, n_tables)

    if storage == "ByteStorage" and shifter == "DefaultShifter":
        return dBG_ByteStorage_DefaultShifter(K, starting_size, n_tables)

    raise TypeError("Invalid Storage or Shifter type.")


def get_dbg_type(str storage='BitStorage',
                 str shifter='DefaultShifter'):
    if storage == "BitStorage" and shifter == "DefaultShifter":
        return dBG_BitStorage_DefaultShifter

    if storage == "NibbleStorage" and shifter == "DefaultShifter":
        return dBG_NibbleStorage_DefaultShifter

    if storage == "ByteStorage" and shifter == "DefaultShifter":
        return dBG_ByteStorage_DefaultShifter

    raise TypeError("Invalid Storage or Shifter type: ({0},{1})".format(storage, shifter))


