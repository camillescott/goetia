from cython.operator cimport dereference as deref

from libc.stdint cimport uint64_t
from libcpp.string cimport string
from libcpp.vector cimport vector

from khmer._oxli.utils import get_n_primes_near_x, is_str, is_num
from boink.utils cimport _bstring, _ustring

{% for Storage_t in Storage_types %}
{% for Shifter_t in Shifter_types %}
{% set tparams %}{{Storage_t}},{{Shifter_t}}{% endset %}
{% set suffix %}{{Storage_t}}_{{Shifter_t}}{% endset %}

cdef class dBG_{{suffix}}:

    def __cinit__(self, int K, uint64_t starting_size, int n_tables):
        cdef vector[uint64_t] primes
        if type(self) is dBG_{{suffix}}:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self._this = make_shared[_dBG[{{tparams}}]](K, primes)

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
        cdef shared_ptr[_KmerIterator[{{Shifter_t}}]] kmer_iter = \
                deref(self._this).get_hash_iter(_sequence)

        cdef hash_t h
        while(not deref(kmer_iter).done()):
            h = deref(kmer_iter).next()
            yield h

    def add_sequence(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef list hits = deref(self._this).add_sequence(_sequence)
        return hits

    # compatibility with oxli API
    consume = add_sequence

    def get_counts(self, str sequence):
        cdef bytes _sequence = _bstring(sequence)
        cdef list counts = deref(self._this).get_counts(_sequence)
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

    def save(self, file_name):
        deref(self._this).save(_bstring(file_name))

    @classmethod
    def load(cls, file_name):
        cdef dBG_{{suffix}} obj = cls(1, 1, 1)
        deref(obj._this).load(_bstring(file_name))
        return obj

{% endfor %}
{% endfor %}


cdef object _make_dbg(int K, uint64_t starting_size, int n_tables,
                      str storage='BitStorage',
                      str shifter='RollingHashShifter'):
    {% for Storage_t in Storage_types %}
    {% set outer_first = loop.first %}
    {% for Shifter_t in Shifter_types %}
    {% set conditional %}storage == "{{Storage_t}}" and shifter == "{{Shifter_t}}"{% endset %}
    {% set suffix %}{{Storage_t}}_{{Shifter_t}}{% endset %}
    {% if outer_first and loop.first %}
    if {{conditional}}:
    {% else %}
    elif {{conditional}}:
    {% endif %}
        return dBG_{{suffix}}(K, starting_size, n_tables)
    {% endfor %}
    {% endfor %}
    else:
        raise TypeError("Invalid Storage or Shifter type.")



