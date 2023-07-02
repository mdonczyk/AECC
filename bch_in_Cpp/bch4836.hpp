#include <atomic>
#include <bit>
#include <bitset>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <memory>
#include <mutex>
#include <random>
#include <set>
#include <thread>
#include <time.h>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#define GFB 63 // Galois Field Barrier = 2**m - 1 = 2**6 - 1
#define n 48
#define k 36
#define t 2
#define HEADER_BYTES 30

// global std::atomic counters:
std::atomic<int> g_success_count{0};
std::atomic<int> g_failure_count{0};
std::atomic<int> g_introduced_errors_count{0};
std::atomic<int> g_big_errors_count{0};
std::atomic<int> g_uncaught_errors_count{0};


// verbose logs and sequential threads flag (off by default):
bool verbose_flag = false;

// stylistic:
const std::string LINE(n/2 - 4, '*');
const std::string DASH_LINE(n/2, '-');

enum Status {
    SUCCESS = 0,
    FAIL = -1
};

class BCH_code_short_t2; // forward declaration

namespace bch {
    // variables:
    int m;
    int primitive_polynomial;
    int alpha_to[GFB], index_of[GFB];
    std::bitset <n> p;
    std::bitset <n> generator_polynomial_bitset;
    std::vector <int> zeros, g, errpos;
    std::vector <std::vector <int>> zeros_cosets;
    std::string filename;
    int error_probability;
    std::mutex Mutex;
    // std::vectors:
    std::vector <std::shared_ptr<BCH_code_short_t2>> BCH_objects;
    std::vector <std::bitset <k>> vector_of_message_polynomials;
    std::vector <unsigned char> recovered_charstream;
    std::vector <unsigned char> modified_charstream;
    // functions:
    /**
     * Read the primitive polynomial of degree 6 from binary form
    */
    void read_p();
    /**
     * Generate a Galois Field for GF(2**m-1) where m = 6
    */
    void generate_gf();
    /**
     * Compute the generator polynomial
    */
    void gen_poly();
    inline int MSB(const std::bitset <n> &Polynomial);
    void verbose_polynomial(const std::bitset <n> &Polynomial);
    int multiply_int_polynomials(int Mulitplicand, int Multiplicator);
    template <size_t N>
    void reverse_bitset(std::bitset <N> &Polynomial, int Shift);
    /**
     * Initialize the primitive polynomial and generate Galois Field and 
     * generator polynomial for later use in encoding and decoding
    */
    void initialize_BCH() {
        read_p();
        generate_gf();
        gen_poly();
    }
};

class BCH_code_short_t2 {
	public:
        explicit BCH_code_short_t2 (const std::bitset <n> &data) : Data(data) {
        }
        /**
         * Calculate redundant bits and encode message into a Codeword polynomial
         */
        void encode_bch();
        /**
         * Introduce errors to a Codeword based on a given probability and save it 
         * as a Received_Codeword
        */
        void introduce_errors();
        /**
         * Print original Codeword and Received_Codeword in binary form and count number of
         * all errors in Received_Codewords and also all errors over t in Received_Codewords
         */
        void print_original_codeword_and_received_codeword();
        /**
         * Calculate syndromes and use Berlekamp-Massey algorith to decode the Received_Codeword
         * polynomial
         * @returns Decoding Status flag
         */
        Status decode_bch();
        /**
         * Print original Data and Decoded_Data in binary form and count number of
         * all uncaught decoding errors
         */
        void print_original_message_and_decoded_message();
        // variables:
        std::bitset <n> Data;
        std::bitset <n> Codeword;
        std::bitset <n> Received_Codeword;
        std::bitset <k> Received_Data;
        std::bitset <k> Decoded_Data;

    private:
        std::vector <int> calculate_syndromes(const std::bitset<n> &Received_Codeword, bool &errors_in_codeword);
        std::pair<std::bitset <n>, std::bitset <n>> divide_bitset_polynomials(const std::bitset <n> &Dividend, const std::bitset <n> &Divisor);
};