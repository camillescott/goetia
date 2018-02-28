cimport cython

from libc.stdint cimport uint8_t, uint16_t, uint64_t

from libcpp cimport bool
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

from boink.hashing cimport *

cdef extern from "oxli/storage.hh":
    # Need these for the Storage template parameter;
    # they don't need methods
    cdef cppclass Storage "oxli::Storage":
        pass
    cdef cppclass BitStorage "oxli::BitStorage" (Storage):
        pass
    cdef cppclass NibbleStorage "oxli::NibbleStorage" (Storage):
        pass
    cdef cppclass QFStorage "oxli::QFStorage" (Storage):
        pass
    cdef cppclass ByteStorage "oxli::ByteStorage" (Storage):
        pass


cdef extern from "boink/dbg.hh" namespace "boink":
    ctypedef pair[bool, bool] bit_pair_t
    ctypedef vector[bit_pair_t] bit_pair_vector_t
    ctypedef uint8_t count_t
    ctypedef pair[count_t, count_t] full_count_t

    cdef cppclass _dBG "boink::dBG" [StorageType, HashShifter] (_KmerClient):
        _dBG(uint16_t, vector[uint64_t])

        hash_t hash(string&)
        bool add(hash_t)
        bool add(string&)

        count_t get(string&)
        count_t get(hash_t)

        vector[bool] add_sequence(string&)

        uint64_t n_unique()
        uint64_t n_occupied()

        uint8_t ** get_raw()

        void save(string)
        void load(string)
        void reset()

    ctypedef _dBG[BitStorage, DefaultShifter] DefaultDBG


include "dbg_types.pxd.pxi"
