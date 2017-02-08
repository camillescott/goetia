from khmer._oxli._oxli cimport CpHashgraph

from khmer._oxli.utils cimport _get_n_primes_near_x


cdef class ProbabilisticGraph:

    def __cinit__(self, int k, int starting_size, int n_tables):
        self.table_sizes = _get_n_primes_near_x(n_tables, starting_size)
        self.n_tables = n_tables
        self.K = k

    cdef 
