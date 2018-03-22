from cython.operator cimport dereference as deref

from libcpp.memory cimport shared_ptr, make_shared
from libcpp.string cimport string

from boink.dbg cimport *
from boink import dbg

include "assembly.pyx.pxi"

def make_assembler(dBG_Base graph):
    return _make_assembler(graph)


def test_cdefs():
    cdef shared_ptr[DefaultDBG] G = make_shared[DefaultDBG](21, 10000, 4)
    cdef shared_ptr[_AssemblerMixin[DefaultDBG]] asm = make_shared[_AssemblerMixin[DefaultDBG]](G)
    
    cdef string S = b"TACATGCTTACTCACAGCACCCCTTCTAATGCGTAACCAGGCGTGGAATT"
    deref(G).add_sequence(S)
    cdef Path path

    deref(asm).assemble(S.substr(0,21), path)
    print(deref(asm).to_string(path))


def test_wrapper():
    G = dbg.make_dBG(21, 10000, 4)
    asm = make_assembler(G)

    S = "TACATGCTTACTCACAGCACCCCTTCTAATGCGTAACCAGGCGTGGAATT"
    G.add_sequence(S)

    print(asm.assemble(S[:21]))
