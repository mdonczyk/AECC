#include <algorithm>
#include <bit>
#include <bitset>
#include <chrono>
#include <cmath>
#include <cstring>
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
#define t 3
#define HEADER_BYTES 30

// global counters:
int success = 0;
int failure = 0;
int uncaught_errors = 0;
int big_errors = 0;
int total_errors = 0;

// verbose logs and sequential threads flag (off by default):
bool verbose = false;

// stylistic:
const string LINE(n/2 - 4, '*');
const string DASH_LINE(n/2, '-');

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
    int error_probability;
    mutex Mutex;
    // vectors:
    vector <unique_ptr<BCH_code_long_t3>> BCH_objects;
    vector <bitset <k>> vector_of_bitsets;
    vector <bitset <k>> vector_of_modified_bitsets;
	vector <bitset <k>> vector_of_recovered_bitsets;
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
        BCH_code_long_t3 (const bitset <n> &data, const int error_probability) : Data(data) {
            Codeword = encode_bch(Data);
            Received_Codeword = introduce_errors(Codeword, error_probability);
            print_original_codeword_and_received_codeword(Codeword, Received_Codeword);
            Received_Data = bitset <k>(Received_Codeword.to_string().substr(0, k));
            Decoded_Data = decode_bch(Received_Codeword);
            // check the decoding status
            if (Decoded_Data.second) {
                print_original_message_and_decoded_message(Data, Decoded_Data.first);
                success++;
            } else {
                failure++;
            }
        }
        /**
         * Calculate redundant bits and encode message into a Codeword polynomial
         * @returns Encoded Data as a Codeword
         */
        bitset <n> encode_bch(const bitset<n> &Data);
        /**
         * Introduce errors to a Codeword based on a given probability and save it 
         * as a Received_Codeword
         * @returns Received_Codeword
        */
        bitset <n> introduce_errors(const bitset <n> &Codeword, const int error_probability);
        /**
         * Print original Codeword and Received_Codeword in binary form and count number
         * all errors in Received_Codewords and also all errors over t in Received_Codewords
         */
        void print_original_codeword_and_received_codeword(const bitset <n> &Codeword, const bitset <n> &Received_Codeword);
        /**
         * Calculate syndromes and use Berlekamp-Massey algorith to decode the Received_Codeword
         * polynomial
         * @returns Decoded message and success flag
         */
        pair<bitset<k>, bool> decode_bch(const bitset <n> &Received_Codeword);
        /**
         * Print original Data and Decoded_Data in binary form and count number
         * all uncaught decoding errors
         */
        void print_original_message_and_decoded_message(const bitset <n> &Data, const bitset <k> &Decoded_Data);
        // variables:
        bitset <n> Data;
        bitset <n> Codeword;
        bitset <n> Received_Codeword;
        bitset <k> Received_Data;
        pair<bitset <k>, bool> Decoded_Data;

    private:
        vector <int> calculate_syndromes(const bitset<n> &Received_Codeword, bool &errors_in_codeword);
        pair<bitset <n>, bitset <n>> divide_bitset_polynomials(const bitset <n> &Dividend, const bitset <n> &Divisor);
};