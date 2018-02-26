#include "boink/hashing.hh"

#include <iostream>

using namespace boink;

namespace boink {

template<class Derived, const std::string& Alphabet>
const std::string HashShifter<Derived, Alphabet>::symbols = Alphabet;

};

int main() {
    std::string S = "ATAT";
    RollingHashShifter<> shifter(S, 4);

    std::cout << shifter.get() << std::endl;
    std::cout << "+++ Standalone hash() +++" << std::endl;
    std::cout << shifter.hash(S) << std::endl;
    std::cout << shifter.hash(S) << std::endl;
    std::cout << "+++ Shifting +++" << std::endl;
    std::cout << shifter.symbol_deque.size() << std::endl;
    std::cout << "cursor: " << shifter.get_cursor() << std::endl;
    std::cout << "shifted hash: " << shifter.shift_left('C') << std::endl;
    std::cout << "shifted cursor: " << shifter.get_cursor() << std::endl;
    std::cout << "standalone shifted: " << shifter.hash(std::string("CATA")) << std::endl;

    std::cout << "shift right G: " << shifter.shift_right('G') << std::endl;
    std::cout << "standalone shifted: " << shifter.hash(std::string("ATAG")) << std::endl;

}
