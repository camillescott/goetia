from boink import libboink
from boink.storage import types

# Because cDBG::Graph is usually first accessed via StreamingCompactor::Compactor, its
# pythonizor seems to be fired after its first instantation already exists. This makes sure
# it's been loaded once already.
for storage_t, args in types:
    _ = libboink.dBG[storage_t, libboink.hashing.RollingHashShifter].build(5, *args)
    _ = libboink.cdbg.cDBG[type(_)].Graph


