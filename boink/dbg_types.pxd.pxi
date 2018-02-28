from libcpp.memory cimport shared_ptr, make_shared


cdef class dBG_BitStorage_DefaultShifter:
    cdef shared_ptr[_dBG[BitStorage,DefaultShifter]] _this
    cdef hash_t _handle_kmer(self, object)

cdef class dBG_NibbleStorage_DefaultShifter:
    cdef shared_ptr[_dBG[NibbleStorage,DefaultShifter]] _this
    cdef hash_t _handle_kmer(self, object)

cdef class dBG_ByteStorage_DefaultShifter:
    cdef shared_ptr[_dBG[ByteStorage,DefaultShifter]] _this
    cdef hash_t _handle_kmer(self, object)

cdef object _make_dbg(int K, uint64_t starting_size, int n_tables, 
                     str storage=*, str shifter=*)