# boink/hashing.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libc.stdint cimport uint8_t, uint16_t, uint64_t

from libcpp cimport bool
from libcpp.deque cimport deque
from libcpp.string cimport string
from libcpp.utility cimport pair
from libcpp.vector cimport vector


cdef extern from "boink/hashing/alphabets.hh" namespace "boink::hashing" nogil:
    # This is a hack to trick Cython into taking static const
    # std::string objects as a template non-type parameters (by default
    # Cython only accepts regular template type parameters)

    cdef cppclass _DNA_SIMPLE   "boink::hashing::DNA_SIMPLE":
        pass
    cdef cppclass _DNAN_SIMPLE  "boink::hashing::DNAN_SIMPLE":
        pass
    cdef cppclass _RNA_SIMPLE   "boink::hashing::RNA_SIMPLE":
        pass
    cdef cppclass _RNAN_SIMPLE  "boink::hashing::RNAN_SIMPLE":
        pass
    cdef cppclass _IUPAC_NUCL   "boink::hashing::IUPAC_NUCL":
        pass
    cdef cppclass _IUPAC_AA     "boink::hashing::IUPAC_AA":
        pass


cdef extern from "boink/hashing/hashing_types.hh" namespace "boink::hashing" nogil:
    ctypedef uint64_t hash_t
    ctypedef pair[hash_t, hash_t] full_hash_t

    cdef struct shift_t:
        shift_t(hash_t, char)
        hash_t hash
        char symbol

    cdef struct kmer_t:
        kmer_t(const hash_t, const string)
        const hash_t hash
        const string kmer

cdef extern from "boink/hashing/hashshifter.hh" namespace "boink::hashing" nogil:

    cdef cppclass _KmerClient "boink::KmerClient":
        _KmerClient(uint16_t)
        uint16_t K()

    cdef cppclass _HashShifter "boink::hashing::HashShifter" [D,A] (_KmerClient):
        hash_t set_cursor(string&)
        string get_cursor()
        void get_cursor(deque[char]&)

        bool is_valid(const char)
        bool is_valid(const string&)

        hash_t get()
        hash_t hash(string&)

        vector[shift_t] gather_left()
        vector[shift_t] gather_right()

        hash_t shift_left(const char)
        hash_t shift_right(const char)

cdef extern from "boink/hashing/rollinghashshifter.hh" namespace "boink::hashing" nogil:
    cdef cppclass _RollingHashShifter "boink::hashing::RollingHashShifter" [A] (_HashShifter[_RollingHashShifter[A], A]):
        _RollingHashShifter(string&, uint16_t)
        _RollingHashShifter(uint16_t)

    ctypedef _RollingHashShifter[_DNA_SIMPLE] _DefaultShifter "boink::hashing::DefaultShifter"

cdef extern from "boink/hashing/kmeriterator.hh" namespace "boink::hashing" nogil:
    cdef cppclass _KmerIterator "boink::hashing::KmerIterator" [S] (_KmerClient):
        _KmerIterator(const string, uint16_t)

        hash_t first()
        hash_t next()
        bool done()

        unsigned int get_start_pos()
        unsigned int get_end_pos()


cdef void _test()
