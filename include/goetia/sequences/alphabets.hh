/**
 * (c) Camille Scott, 2019
 * File   : alphabets.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 07.01.2020
 */


#ifndef GOETIA_ALPHABETS_HH
#define GOETIA_ALPHABETS_HH

#include <sstream>
#include <string>
#include <string_view>

#include "goetia/sequences/exceptions.hh"

#define rc_tbl \
  "                                                                "\
  /*ABCDEFGHIJKLMNOPQRSTUVWXYZ      abcdefghijklmnopqrstuvwxyz    */\
  " TVGH FCD  M KN   YSAABW R       TVGH FCD  M KN   YSAABW R"
  //" TVGH FCD  M KA   YSAABWARA      TVGH FCD  M KA   YSAABWARA"

namespace goetia {

template <class Derived>
struct Alphabet {

    static constexpr auto SYMBOLS = std::string_view("NNNN");
    static constexpr auto COMPLEMENTS = std::string_view("NNNN");

    static const size_t size() {
        return SYMBOLS.size();
    }

    static const char validate(const char c) {
        return Derived::_validate(c);
    }

    static void validate(char * sequence, const size_t length) {
        for (size_t i = 0; i < length; ++i) {
            const char validated = validate(sequence[i]);
            if (validated == '\0') {
                std::ostringstream os;
                os << "Alphabet: Invalid symbol '"
                   << sequence[i] << "' in sequence "
                   << std::string(sequence, length)
                   << " (alphabet=" << SYMBOLS << ").";

                throw InvalidCharacterException(os.str().c_str());
            } else {
                sequence[i] = validated;
            }
        }
    }

    static void validate(const char * sequence, const size_t length) {
        for (size_t i = 0; i < length; ++i) {
            const char validated = validate(sequence[i]);
            if (validated == '\0') {
                std::ostringstream os;
                os << "Alphabet: Invalid symbol '"
                   << sequence[i] << "' in sequence "
                   << std::string(sequence, length)
                   << " (alphabet=" << SYMBOLS << ").";

                throw InvalidCharacterException(os.str().c_str());
            } 
        }
    }
    
    static const char complement(const char c) {
        return Derived::_complement(c);
    }

    static std::string reverse_complement(const std::string& sequence) {
        return Derived::_reverse_complement(sequence);
    }

private:

    friend Derived;

};


struct ANY : public Alphabet<ANY> {

    static constexpr auto SYMBOLS = std::string_view("*");
    static constexpr auto COMPLEMENTS = std::string_view("*");

    static const char _validate(const char c) {
        return toupper(c);
    }

    static const char _complement(const char c) {
        return c;
    }

    static std::string _reverse_complement(const std::string& sequence) {
        return std::string(sequence);
    }
};


struct DNA_SIMPLE : public Alphabet<DNA_SIMPLE> {

    static constexpr auto SYMBOLS = std::string_view("ACGT");
    static constexpr auto COMPLEMENTS = std::string_view("TGCA");

    static const char _validate(const char c) {
        switch(c) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
                return c;
            case 'a':
                return 'A';
            case 'c':
                return 'C';
            case 'g':
                return 'G';
            case 't':
                return 'T';
            default:
                return '\0';
        }
    }

    static const char _complement(const char c) {
        switch(c) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            default:
                return '\0';
        }
    }

    static std::string _reverse_complement(const std::string& sequence) {
        std::string out = sequence;

        auto from = out.begin();
        auto to = out.end();

        char c;
        for (to--; from <= to; from++, to--) {
            c = rc_tbl[(int)*from];
            *from = rc_tbl[(int)*to];
            *to = c;
        }
        return out;
    }
};

struct DNAN_SIMPLE : public Alphabet<DNAN_SIMPLE> {

    static constexpr auto SYMBOLS = std::string_view("ACGTN");
    static constexpr auto COMPLEMENTS = std::string_view("TGCAN");

    static const char _validate(const char c) {
        switch(c) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'N':
                return c;
            case 'a':
                return 'A';
            case 'c':
                return 'C';
            case 'g':
                return 'G';
            case 't':
                return 'T';
            case 'n':
                return 'N';
            default:
                return '\0';
        }
    }

    static const char _complement(const char c) {
        switch(c) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            case 'N':
                return 'N';
            default:
                return '\0';
        }
    }

    static std::string _reverse_complement(const std::string& sequence) {
        std::string out = sequence;

        auto from = out.begin();
        auto to = out.end();

        char c;
        for (to--; from <= to; from++, to--) {
            c = _complement(*from);
            *from = _complement(*to);
            *to = c;
        }
        return out;
    }
};


struct IUPAC_NUCL : public Alphabet<IUPAC_NUCL> {

    // ref: http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html#iupacdegeneracies
    static constexpr auto SYMBOLS = "ATUGCYRSWKMBDHVN";
    static constexpr auto COMPLEMENTS = "TAACGRYSWMKVHDBN";
    
    static const char _validate(const char c) {
        switch(c) {
            case 'A':
            case 'T':
            case 'U':
            case 'G':
            case 'C':
            case 'Y':
            case 'R':
            case 'S':
            case 'W':
            case 'K':
            case 'M':
            case 'B':
            case 'D':
            case 'H':
            case 'V':
            case 'N':
                return c;
            case 'a':
                return 'A';
            case 't':
                return 'T';
            case 'u':
                return 'U';
            case 'g':
                return 'G';
            case 'c':
                return 'C';
            case 'y':
                return 'Y';
            case 'r':
                return 'R';
            case 's':
                return 'S';
            case 'w':
                return 'W';
            case 'k':
                return 'K';
            case 'm':
                return 'M';
            case 'b':
                return 'B';
            case 'd':
                return 'D';
            case 'h':
                return 'H';
            case 'v':
                return 'V';
            case 'n':
                return 'N';
            default:
                return '\0';
        }
    }

    static const char _complement(const char c) {
        switch(c) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'U':
                return 'A';
            case 'G':
                return 'C';
            case 'C':
                return 'G';
            case 'Y':
                return 'R';
            case 'R':
                return 'Y';
            case 'S':
                return 'S';
            case 'W':
                return 'W';
            case 'K':
                return 'M';
            case 'M':
                return 'K';
            case 'B':
                return 'V';
            case 'D':
                return 'H';
            case 'H':
                return 'D';
            case 'V':
                return 'B';
            case 'N':
                return 'N';
            default:
                return '\0';
        }
    }

    static std::string _reverse_complement(const std::string& sequence) {
        std::string out = sequence;

        auto from = out.begin();
        auto to = out.end();

        char c;
        for (to--; from <= to; from++, to--) {
            c = _complement(*from);
            *from = _complement(*to);
            *to = c;
        }
        return out;
    }
};

}

#endif
