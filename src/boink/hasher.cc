#include "boink/hasher.hh"

using namespace boink;

namespace boink {

template<class Derived, const std::string& Alphabet>
const std::string HashShifter<Derived, Alphabet>::symbols = Alphabet;

};
