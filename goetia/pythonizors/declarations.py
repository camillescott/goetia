from goetia import libgoetia
from goetia.storage import _types as storage_types
from goetia.hashing import types as hashing_types

from cppyy.gbl import std


# Because cDBG::Graph is usually first accessed via StreamingCompactor::Compactor, its
# pythonizor seems to be fired after its first instantation already exists. This makes sure
# it's been loaded once already.
for storage_t, args in storage_types:
    for hasher_t in hashing_types:
        if 'UnikmerLemirePolicy' in hasher_t.__name__:
            hasher = hasher_t.build(21, 7)
        else:
            hasher = hasher_t(21)
        storage = storage_t.build(*args)
        dbg = libgoetia.dBG[storage_t, hasher_t].build(storage, hasher)
        _ = libgoetia.cdbg.cDBG[type(dbg)].CompactNode
        _ = libgoetia.cdbg.cDBG[type(dbg)].UnitigNode
        _ = libgoetia.cdbg.cDBG[type(dbg)].DecisionNode
        _ = libgoetia.cdbg.cDBG[type(dbg)].Graph

