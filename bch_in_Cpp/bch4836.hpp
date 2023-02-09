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
#define n 48
#define t 2
#define k 36
#define HEADER_BYTES 30
typedef unsigned long long ULL;

namespace BCH {
    // static variables:
    static int primitive_polynomial;
    static int alpha_to[GF], index_of[GF];
    static bitset <n> p;
    static bitset <n> generator_polynomial_bitset;
    static vector <int> zeros, g, errpos;
    static vector <vector <int>> zeros_deluxe;
    // functions:
    void read_p();
    void generate_gf();
    void gen_poly();
    template <size_t N>
    int MSB(const bitset <N> &Polynomial);
    void verbose_polynomial(const bitset <n> &Polynomial);
    uint multiply_uint_polynomials(uint Mulitplicand, uint Multiplicator);
    template <size_t N>
    void reverse_bitset(bitset <N> &Polynomial, int Shift);
    void initialize_BCH() {
        read_p();
        generate_gf();
        gen_poly();
    }
};

class BCH_code_short_t2 {
	public:
        bitset <n> encode_bch(const bitset <n> &Data);
        void print_codeword_and_received_codeword(const bitset <n> &Codeword, const bitset <n> &Received_Codeword);
        void print_message_and_decoded_message(const bitset <n> &Data, const bitset <k> &Decoded_Data);
        bitset <n> introduce_errors(const bitset <n> &Codeword, const int &Probability);
        pair<bitset <k>, bool> decode_bch(const bitset <n> &Received_Codeword);
        // variables
        bitset<n> Data;
        bitset<n> Codeword;
        bitset<n> Received_Codeword;
        bitset <k> Received_Data;
        pair<bitset <k>, bool> Decoded_Data;

    private:
        vector <int> calculate_syndromes(const bitset <n> &Received_Codeword, bool &Syn_error);;
        bitset <n> multiply_bitset_polynomials(const bitset <n> &Mulitplicand, const bitset <n> &Multiplicator);
        pair<bitset <n>, bitset <n>> divide_bitset_polynomials(const bitset <n> &Dividend, const bitset <n> &Divisor);
};
