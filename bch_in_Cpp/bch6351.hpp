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
#include <random>
#include <set>
#include <time.h>
#include <vector>

using namespace std;

#define GF 63 // Galois Field --> 2**m - 1 = 2**6 - 1
#define n 63
#define t 2
#define k 51
#define HEADER_BYTES 30
typedef unsigned long long ULL;

// global counters:
int success = 0;
int failure = 0;
int uncaught_errors = 0;
int big_errors = 0;
int total_errors = 0;
bool verbose = false;

// stylistic:
string line(n/2, '*');
string dash_line(n/2, '-');

namespace BCH {
    // static variables:
    static int primitive_polynomial;
    static int alpha_to[GF], index_of[GF];
    static bitset <n> p;
    static bitset <n> generator_polynomial_bitset;
    static vector <int> zeros, g, errpos;
    static vector <vector <int>> zeros_deluxe;
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
    template <size_t N>
    int MSB(const bitset <N> &Polynomial);
    void verbose_polynomial(const bitset <n> &Polynomial);
    uint multiply_uint_polynomials(uint Mulitplicand, uint Multiplicator);
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

class BCH_code_long_t2 {
	public:
        BCH_code_long_t2 (bitset<n> data, int error_probability) {
            Data = data;
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
        vector <int> calculate_syndromes(bool &Syn_error);;
        bitset <n> multiply_bitset_polynomials(const bitset <n> &Mulitplicand, const bitset <n> &Multiplicator);
        pair<bitset <n>, bitset <n>> divide_bitset_polynomials(const bitset <n> &Dividend, const bitset <n> &Divisor);
};
