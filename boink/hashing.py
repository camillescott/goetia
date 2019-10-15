from boink import libboink


_types = [libboink.hashing.RollingHashShifter,
          libboink.hashing.UKHS.LazyShifter]
types = {hashing_t.__name__ : hashing_t for hashing_t in _types}


RollingHashShifter = libboink.hashing.RollingHashShifter

UKHS               = libboink.hashing.UKHS
Unikmer            = UKHS.Unikmer
BinnedKmer         = UKHS.BinnedKmer
UnikmerMap         = UKHS.Map
UnikmerShifter     = UKHS.LazyShifter
