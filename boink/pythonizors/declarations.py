from boink import libboink
from boink.storage import _types as storage_types
from boink.hashing import types as hashing_types

from cppyy.gbl import std


# Because cDBG::Graph is usually first accessed via StreamingCompactor::Compactor, its
# pythonizor seems to be fired after its first instantation already exists. This makes sure
# it's been loaded once already.
for storage_t, args in storage_types:
    for hasher_t in hashing_types:
        if hasher_t.__name__.startswith('UnikmerShifter'):
            hasher = hasher_t.build(21, 7)
        else:
            hasher = hasher_t(21)
        storage = storage_t.build(*args)
        #hasher.set_cursor('A' * 21)
        dbg = libboink.dBG[storage_t, hasher_t].build(storage, hasher)
        _ = libboink.cdbg.cDBG[type(dbg)].CompactNode
        _ = libboink.cdbg.cDBG[type(dbg)].UnitigNode
        _ = libboink.cdbg.cDBG[type(dbg)].DecisionNode
        _ = libboink.cdbg.cDBG[type(dbg)].Graph

