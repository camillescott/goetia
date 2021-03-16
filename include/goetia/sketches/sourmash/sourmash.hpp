#ifndef KMER_MIN_HASH_HH
#define KMER_MIN_HASH_HH

#include <algorithm>
#include <exception>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <string>


namespace sourmash {

    extern "C" {
      #include "goetia/sketches/sourmash/sourmash.h"
    }

    inline uint64_t _hash_murmur(const std::string& kmer, const uint32_t seed) {
        return hash_murmur(kmer.c_str(), seed);
    }

    typedef uint64_t HashIntoType;

    class sourmash_exception : public std::exception {
      public:
        explicit sourmash_exception(const std::string& msg = "Generic sourmash.rs exception") : _msg(msg) {
        }

        virtual ~sourmash_exception() throw() {
        }
        virtual const char* what() const throw() {
            return _msg.c_str();
        }

      protected:
        const std::string _msg;
    };

    inline void process_errors() {
        auto err_code = sourmash_err_get_last_code();

        switch (err_code) {

            case SOURMASH_ERROR_CODE_NO_ERROR:
                break;
            case SOURMASH_ERROR_CODE_MISMATCH_K_SIZES:
                throw sourmash_exception("sourmash.rs: different ksizes cannot be compared");
                break;
            case SOURMASH_ERROR_CODE_MISMATCH_DNA_PROT:
                throw sourmash_exception("sourmash.rs: DNA/prot minhashes cannot be compared");
                break;
            case SOURMASH_ERROR_CODE_MISMATCH_SCALED:
                throw sourmash_exception("sourmash.rs: mismatch in scaled param; comparison fail");
                break;
            case SOURMASH_ERROR_CODE_MISMATCH_SEED:
                throw sourmash_exception("sourmash.rs: mismatch in seed; comparison fail");
                break;
            case SOURMASH_ERROR_CODE_INVALID_DNA:
                throw sourmash_exception("sourmash.rs: invalid DNA given to kmerminhash_add_sequence");
            default:
                throw sourmash_exception(std::string("sourmash.rs error: code ") + std::to_string(err_code));
                break;
        }
    }

    class MinHash {
      protected:
        SourmashKmerMinHash* _this;

      public:
        MinHash(unsigned int n, unsigned int k, bool prot, bool dayhoff, bool hp, uint32_t s, HashIntoType mx) {
            _this = kmerminhash_new(n, k, prot, dayhoff, hp, s, mx, false);
        };

        MinHash(SourmashKmerMinHash* _this) : _this(_this) {
        }

        void add_hash(const HashIntoType h) {
            kmerminhash_add_hash(_this, h);
        }

        void remove_hash(const HashIntoType h) {
            kmerminhash_remove_hash(_this, h);
        }

        void add_word(const std::string& word) {
            kmerminhash_add_word(_this, word.c_str());
        }

        void add_sequence(const char* sequence, bool force = false) {
            kmerminhash_add_sequence(_this, sequence, force);
            process_errors();
        }

        void merge(const MinHash& other) {
            kmerminhash_merge(_this, other._this);
            process_errors();
        }

        unsigned int count_common(const MinHash& other, bool downsample = false) {
            auto v = kmerminhash_count_common(_this, other._this, downsample);
            process_errors();
            return v;
        }

        size_t size() {
            return kmerminhash_get_mins_size(_this);
        }

        uint32_t num() {
            return kmerminhash_num(_this);
        }

        uint64_t seed() {
            return kmerminhash_seed(_this);
        }

        bool track_abundance() {
            return kmerminhash_track_abundance(_this);
        }

        bool is_protein() {
            return kmerminhash_is_protein(_this);
        }

        bool dayhoff() {
            return kmerminhash_dayhoff(_this);
        }

        bool hp() {
            return kmerminhash_hp(_this);
        }

        uint32_t ksize() {
            return kmerminhash_ksize(_this);
        }

        uint64_t max_hash() {
            return kmerminhash_max_hash(_this);
        }

        char aa_to_dayhoff(char aa) {
            return sourmash_aa_to_dayhoff(aa);
        }

        char translate_codon(const char* codon) {
            return sourmash_translate_codon(codon);
        }

        std::vector<HashIntoType> mins() {
            uintptr_t * size = new uintptr_t;
            auto        ptr = kmerminhash_get_mins(_this, size);
            std::vector<HashIntoType> m(ptr, ptr + *size);
            delete size;
            return m;
        }

        ~MinHash() throw() {
            kmerminhash_free(_this);
        }

        SourmashKmerMinHash* _get_ptr() {
            return _this;
        }
    };

    class MinAbundance : public MinHash {
      public:
        MinAbundance(unsigned int n, unsigned int k, bool prot, bool dayhoff, bool hp, uint32_t seed, HashIntoType mx)
            : MinHash(n, k, prot, dayhoff, hp, seed, mx) {
            kmerminhash_free(_this);
            _this = kmerminhash_new(n, k, prot, dayhoff, hp, seed, mx, true);
        };

        std::vector<HashIntoType> abunds() {
            uintptr_t * size = new uintptr_t;
            auto        ptr = kmerminhash_get_abunds(_this, size);
            std::vector<HashIntoType> m(ptr, ptr + *size);
            delete size;
            return m;
        }

        void set_abundances(std::vector<HashIntoType> mins, std::vector<HashIntoType> abunds, bool clear = false) {
            // TODO: assert mins and abunds are the same size?
            kmerminhash_set_abundances(_this, mins.data(), abunds.data(), mins.size(), clear);
        }

        ~MinAbundance() throw() {
        }
    };

}  // namespace sourmash

#endif  // KMER_MIN_HASH_HH
