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

using namespace std;

#define GFB 63 // Galois Field Barrier = 2**m - 1 = 2**6 - 1
#define n 63
#define k 45
#define t 2
#define HEADER_BYTES 30

// global atomic counters:
atomic<int> g_success_count{0};
atomic<int> g_failure_count{0};
atomic<int> g_introduced_errors_count{0};
atomic<int> g_big_errors_count{0};
atomic<int> g_uncaught_errors_count{0};


// verbose logs and sequential threads flag (off by default):
bool verbose_flag = false;

// stylistic:
const string LINE(n/2 - 4, '*');
const string DASH_LINE(n/2, '-');

enum Status {
    SUCCESS = 0,
    FAIL = -1
};

class BCH_code_long_t3; // forward declaration

namespace BCH {
    // variables:
    int m;
    int primitive_polynomial;
    int alpha_to[GFB], index_of[GFB];
    bitset <n> p;
    bitset <n> generator_polynomial_bitset;
    vector <int> zeros, g, errpos;
    vector <vector <int>> zeros_cosets;
    string filename;
    int error_probability;
    mutex Mutex;
    // vectors:
    vector <shared_ptr<BCH_code_long_t3>> BCH_objects;
    vector <bitset <k>> vector_of_message_polynomials;
    vector <unsigned char> recovered_charstream;
    vector <unsigned char> modified_charstream;
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
    int MSB(const bitset <n> &Polynomial);
    void verbose_polynomial(const bitset <n> &Polynomial);
    int multiply_int_polynomials(int Mulitplicand, int Multiplicator);
    template <size_t N>
    void reverse_bitset(bitset <N> &Polynomial, int Shift);
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

class BCH_code_long_t3 {
	public:
        explicit BCH_code_long_t3 (const bitset <n> &data) : Data(data) {
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
        bitset <n> Data;
        bitset <n> Codeword;
        bitset <n> Received_Codeword;
        bitset <k> Received_Data;
        bitset <k> Decoded_Data;

    private:
        vector <int> calculate_syndromes(const bitset<n> &Received_Codeword, bool &errors_in_codeword);
        pair<bitset <n>, bitset <n>> divide_bitset_polynomials(const bitset <n> &Dividend, const bitset <n> &Divisor);
};