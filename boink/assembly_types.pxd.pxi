from libcpp.memory cimport shared_ptr, make_shared

from boink.dbg cimport *

cdef class Assembler_Base:
    cdef readonly object storage_type
    cdef readonly object shifter_type



cdef class Assembler_BitStorage_DefaultShifter(Assembler_Base):
    cdef shared_ptr[_AssemblerMixin[_dBG[BitStorage,DefaultShifter]]] _this
    cdef shared_ptr[_dBG[BitStorage,DefaultShifter]] _graph
    cdef readonly dBG_BitStorage_DefaultShifter Graph



cdef class Assembler_NibbleStorage_DefaultShifter(Assembler_Base):
    cdef shared_ptr[_AssemblerMixin[_dBG[NibbleStorage,DefaultShifter]]] _this
    cdef shared_ptr[_dBG[NibbleStorage,DefaultShifter]] _graph
    cdef readonly dBG_NibbleStorage_DefaultShifter Graph



cdef class Assembler_ByteStorage_DefaultShifter(Assembler_Base):
    cdef shared_ptr[_AssemblerMixin[_dBG[ByteStorage,DefaultShifter]]] _this
    cdef shared_ptr[_dBG[ByteStorage,DefaultShifter]] _graph
    cdef readonly dBG_ByteStorage_DefaultShifter Graph


cdef object _make_assembler(dBG_Base graph)