from boink import libboink


_types = [libboink.hashing.RollingHashShifter]
types = {hashing_t.__name__ : hashing_t for hashing_t in _types}
