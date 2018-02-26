cimport cython

from libc.stdint cimport uint8_t, uint16_t, uint64_t

from libcpp cimport bool
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

from boink.hashing cimport _KmerClient, hash_t

cdef extern from "boink/dbg.hh":
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
        count_t get(hash_t&)

        vector[bool] add_sequence(string&)

        void save(string)
        void load(string)

    cdef cppclass _DefaultDBG "boink::DefaultDBG" (_dBG):
        pass

