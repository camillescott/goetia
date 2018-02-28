from cython.operator cimport dereference as deref

from libc.stdint cimport uint64_t
from libcpp.string cimport string
from libcpp.vector cimport vector

from khmer._oxli.utils import get_n_primes_near_x, is_str, is_num
from boink.utils cimport _bstring, _ustring


cdef class dBG_BitStorage_DefaultShifter:

    def __cinit__(self, int K, uint64_t starting_size, int n_tables):
        cdef vector[uint64_t] primes
        if type(self) is dBG_BitStorage_DefaultShifter:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self._this = make_shared[_dBG[BitStorage,DefaultShifter]](K, primes)

    cdef hash_t _handle_kmer(self, object kmer):
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

    def add_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        return deref(self._this).add_sequence(_sequence)

    @property
    def n_unique(self):
        return deref(self._this).n_unique()

    def save(self, file_name):
        deref(self._this).save(_bstring(file_name))

    @classmethod
    def load(cls, file_name):
        cdef dBG_BitStorage_DefaultShifter obj = cls(1, 1, 1)
        deref(obj._this).load(_bstring(file_name))
        return obj


cdef class dBG_NibbleStorage_DefaultShifter:

    def __cinit__(self, int K, uint64_t starting_size, int n_tables):
        cdef vector[uint64_t] primes
        if type(self) is dBG_NibbleStorage_DefaultShifter:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self._this = make_shared[_dBG[NibbleStorage,DefaultShifter]](K, primes)

    cdef hash_t _handle_kmer(self, object kmer):
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

    def add_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        return deref(self._this).add_sequence(_sequence)

    @property
    def n_unique(self):
        return deref(self._this).n_unique()

    def save(self, file_name):
        deref(self._this).save(_bstring(file_name))

    @classmethod
    def load(cls, file_name):
        cdef dBG_NibbleStorage_DefaultShifter obj = cls(1, 1, 1)
        deref(obj._this).load(_bstring(file_name))
        return obj


cdef class dBG_ByteStorage_DefaultShifter:

    def __cinit__(self, int K, uint64_t starting_size, int n_tables):
        cdef vector[uint64_t] primes
        if type(self) is dBG_ByteStorage_DefaultShifter:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self._this = make_shared[_dBG[ByteStorage,DefaultShifter]](K, primes)

    cdef hash_t _handle_kmer(self, object kmer):
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

    def add_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        return deref(self._this).add_sequence(_sequence)

    @property
    def n_unique(self):
        return deref(self._this).n_unique()

    def save(self, file_name):
        deref(self._this).save(_bstring(file_name))

    @classmethod
    def load(cls, file_name):
        cdef dBG_ByteStorage_DefaultShifter obj = cls(1, 1, 1)
        deref(obj._this).load(_bstring(file_name))
        return obj



cdef object _make_dbg(int K, uint64_t starting_size, int n_tables,
                      str storage='BitStorage',
                      str shifter='RollingHashShifter'):
    if storage == "BitStorage" and shifter == "DefaultShifter":
        return dBG_BitStorage_DefaultShifter(K, starting_size, n_tables)
    elif storage == "NibbleStorage" and shifter == "DefaultShifter":
        return dBG_NibbleStorage_DefaultShifter(K, starting_size, n_tables)
    elif storage == "ByteStorage" and shifter == "DefaultShifter":
        return dBG_ByteStorage_DefaultShifter(K, starting_size, n_tables)
    else:
        raise TypeError("Invalid Storage or Shifter type.")


