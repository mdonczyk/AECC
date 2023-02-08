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

class BCH_code_long_t2 {
	public:
        BCH_code_long_t2(){
            read_p();		// read primitive polynomial p(x) 
            generate_gf();	// generate the Galois Field GF(2**m) (GF(64))
            gen_poly();		// Compute the generator polynomial g(x) of BCH code
        }
        vector <bool> bytes_to_bits (const vector <unsigned char> &Buffer);
        vector <char> bits_to_bytes (const vector <bool> &Recovered_bits, const int fileSize);
        vector <bitset <k>> bits_to_bitsets (const vector <bool> &Buffer_bits);
        vector <bool> bitset_to_bits (const vector <bitset <k>> &Bector_of_bits);
        bitset <n> generate_data();
        bitset <n> encode_bch(const bitset <n> &Data);
        void print_codeword_and_received_codeword(const bitset <n> &Codeword, const bitset <n> &Received_Codeword);
        void print_message_and_decoded_message(const bitset <n> &Data, const bitset <k> &Decoded_Data);
        bitset <n> introduce_errors(const bitset <n> &Codeword, const int &Probability);
        template <size_t N>
        string print_wihtout_zeros(const bitset <N> &Polynomial, const uint &Not_Zeros);
        pair<bitset <k>, bool> decode_bch(const bitset <n> &Received_Codeword);
        bitset <n> user_input(const bitset <n> &Codeword);
        void cin_clean();

    private:
        //variables:
        int primitive_polynomial = 0;
        int alpha_to[GF] = {0}, index_of[GF] = {0};
        bitset <n> p;
        bitset <n> generator_polynomial_bitset;
		vector <int> zeros, g, errpos;
		vector <vector <int>> zeros_deluxe;
        //functions:
        void read_p();
        void generate_gf();
        void gen_poly();
        template <size_t N>
        int MSB(const bitset <N> &Polynomial);
        vector <int> calculate_syndromes(const bitset <n> &Received_Codeword, bool &Syn_error);
        void verbose_polynomial(const bitset <n> &Polynomial);
        template <size_t N>
        void reverse_bitset(bitset <N> &Polynomial, int Shift);
        uint multiply_uint_polynomials(uint Mulitplicand, uint Multiplicator);
        bitset <n> multiply_bitset_polynomials(const bitset <n> &Mulitplicand, const bitset <n> &Multiplicator);
        pair<bitset <n>, bitset <n>> divide_bitset_polynomials(const bitset <n> &Dividend, const bitset <n> &Divisor);
};
