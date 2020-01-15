/**
 * (c) Camille Scott, 2019
 * File   : alphabets.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 07.01.2020
 */


#ifndef BOINK_ALPHABETS_HH
#define BOINK_ALPHABETS_HH

#include <sstream>
#include <string>

#include "boink/sequences/exceptions.hh"

#define rc_tbl \
  "                                                                "\
  /*ABCDEFGHIJKLMNOPQRSTUVWXYZ      abcdefghijklmnopqrstuvwxyz    */\
  " TVGH FCD  M KN   YSAABW R       TVGH FCD  M KN   YSAABW R"
  //" TVGH FCD  M KA   YSAABWARA      TVGH FCD  M KA   YSAABWARA"

namespace boink {

template <class Derived>
struct Alphabet {

    static const std::string SYMBOLS;
    static const std::string COMPLEMENTS;

    static const size_t size() {
        return SYMBOLS.size();
    }

    static const char validate(const char c) {
        return Derived::_validate(c);
    }

    static void validate(char * sequence, const size_t length) {
        for (size_t i = 0; i < length; ++i) {
            const char validated = sequence[i];
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
            const char validated = sequence[i];
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


struct DNA_SIMPLE : public Alphabet<DNA_SIMPLE> {

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
