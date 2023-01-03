#include <iostream>
#include <algorithm>
#include <vector>
#include <bitset>
#include <time.h>
#include <cmath>
#include <chrono>
#include <bit>
#include <set>

using namespace std;

#define GF 63 // Galois Field --> 2**m - 1 = 2**6 - 1
#define k 51
#define n 63
#define t 2

class BCH_code {
	public:
        BCH_code(){
            auto start = chrono::high_resolution_clock::now();
            read_p();		// read primitive polynomial p(x) 
            generate_gf();	// generate the Galois Field GF(2**m) (GF(64))
            gen_poly();		// Compute the generator polynomial g(x) of BCH code
            auto stop = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
            cout <<"Count: "<< duration.count() << " microseconds" << endl;
        }
        bitset <GF> generate_data();
        bitset <GF> encode_bch(const bitset <GF> &Data);
        void print_codeword_and_received_codeword(const bitset <GF> &Codeword, const bitset <GF> &Received_Codeword);
        void print_message_and_decoded_message(const bitset <GF> &Data, const bitset <GF> &Decoded_Data);
        string print_wihtout_zeros(const bitset <GF> &Polynomial, const uint &Not_Zeros);
        bitset <GF> decode_bch(const bitset <GF> &Received_Codeword);
        bitset <GF> user_input(const bitset <GF> &Codeword);
        void cin_clean();

    private:
        //variables:
        int primitive_polynomial = 0;
        int alpha_poly_from_index[64] = {0}, index_of_alpha_from_poly[64] = {0};
        int c[GF] = {0}, decerror = 0;
        bitset <GF> p;
        bitset <GF> generator_polynomial_bitset;
		vector <int> zeros, g, errpos;
		vector <vector <int>> zeros_deluxe;
        //functions:
        void read_p();
        void generate_gf();
        void gen_poly();
        int MSB(const bitset <GF> &polynomial);
        vector <int> calculate_syndromes(const bitset <GF> &Received_Codeword, bool &syn_error);
        void verbose_polynomial(const bitset <GF> &polynomial);
        uint multiply_uint_polynomials(uint mulitplicand, uint multiplicator);
        bitset <GF> multiply_bitset_polynomials(const bitset <GF> &mulitplicand, const bitset <GF> &multiplicator);
        pair<bitset <GF>, bitset <GF>> divide_bitset_polynomials(const bitset <GF> &dividend, const bitset <GF> &divisor);
};
		/*
		block n n = 63 --> 64-1  Gf(2**6)
		generate GF(2**m) from the irreducible polynomial p(X) in p[0]..p[m]
		lookup tables:  index->polynomial form   power_of_alpha[] contains j=alpha_poly_from_index**i;
		polynomial form -> index form  alpha_poly_from_index[j=alpha_poly_from_index**i] = i alpha_poly_from_index=2 is the
		primitive element of GF(2**m) 
			
							GF(64) : P(x) = x6 + x4 + x3 + x + 1 = 1011011 = 91

		alpha_poly_from_index[4] = 16)	index_of_alpha_from_poly[4] = 2) 	          	  
-----------------------------------------------------------------------------------------------------------------
index of alpha_poly_from_index:		power of alpha_poly_from_index		     |  MINIMAL POLYNOMIAL (addition of all elements in a given cyclotomic set)
a0	= 1					= 000001		= 1		= a63	= a126	= a189	...	 |  x+1 
a1	= a	 				= 000010		= 2		= a64	= a127	= a190	...	 |  x6+x4+x3+x+1 = 1011011
a2	= a2				= 000100		= 4		= a65	= a128	= a191	...	 |  x6+x4+x3+x+1
a3	= a3				= 001000		= 8		= a66	= a129	= a192	...	 |  x6+x5+x4+x2+1 = 1110101
a4	= a4				= 010000		= 16	= a67	= a130	= a193	...	 |  x6+x4+x3+x+1
a5	= a5				= 100000		= 32	= a68	= a131	= a194	...	 |  x6+x+1
a6	= a4+a3+a+1			= 011011		= 27	= a69	= a132	= a195	...	 |  x6+x5+x4+x2+1
a7	= a5+a4+a2+a		= 110110		= 54	= a70	= a133	= a196	...	 |  x6+x3+1
a8	= a5+a4+a2+a+1		= 110111		= 55	= a71	= a134	= a197	...	 |  x6+x4+x3+x+1
a9	= a5+a4+a2+1		= 110101		= 53	= a72	= a135	= a198	...	 |  x3+x+1
a10	= a5+a4+1			= 110001		= 49	= a73	= a136	= a199	...	 |  x6+x+5+x2+x+1
a11	= a5+a4+a3+1		= 111001		= 57	= a74	= a137	= a200	...	 |  x6+x1
a12	= a5+a3+1			= 101001		= 41	= a75	= a138	= a201	...	 |  x6+x5+x4+x2+1
a13	= a3+1				= 001001		= 9		= a76	= a139	= a202	...	 |  x6+x5+x4+x+1
a14	= a4+a				= 010010		= 18	= a77	= a140	= a203	...	 |  x6+x3+1
a15	= a5+a2				= 100100		= 36	= a78	= a141	= a204	...	 |  x6+x4+x2+x+1
a16	= a4+a+1			= 010011		= 19	= a79	= a142	= a205	...	 |  x6+x4+x3+x+1
a17	= a5+a2+a			= 100110		= 38	= a80	= a143	= a206	...	 |  x6+x+1
a18	= a4+a2+a+1			= 010111		= 23	= a81	= a144	= a207	...	 |  x3+x+1
a19	= a5+a3+a2+a		= 101110		= 46	= a82	= a145	= a208	...	 |  x6+x5+x4+x+1
a20	= a2+a+1			= 000111		= 7		= a83	= a146	= a209	...	 |  x6+x+1
a21	= a3+a2+a			= 001110		= 14	= a84	= a147	= a210	...	 |  x2+x+1
a22	= a4+a3+a2			= 011100		= 28	= a85	= a148	= a211	...	 |  x6+x5+x2+x+1
a23	= a5+a4+a3			= 111000		= 56	= a86	= a149	= a212	...	 |  x6+x5+1
a24	= a5+a3+a+1			= 101011		= 43	= a87	= a150	= a213	...	 |  x6+x5+x4+x2+1
a25	= a3+a2+1			= 001101		= 13	= a88	= a151	= a214	...	 |  x6+x5+x2+x+1
a26	= a4+a3+a			= 011010		= 26	= a89	= a152	= a215	...	 |  x6+x5+x4+x+1
a27	= a5+a4+a2			= 110100		= 52	= a90	= a153	= a216	...	 |  x3+x2+1
a28	= a5+a4+a+1			= 110011		= 51	= a91	= a154	= a217	...	 |  x6+x3+1
a29	= a5+a4+a3+a2+1		= 111101		= 61	= a92	= a155	= a218	...	 |  x6+x5+1
a30	= a5+1				= 100001		= 33	= a93	= a156	= a219	...	 |  x6+x4+x2+x+1
a31	= a4+a3+1			= 011001		= 25	= a94	= a157	= a220	...	 |  x6+x5+x3+x2+1
a32	= a5+a4+a			= 110010		= 50	= a95	= a158	= a221	...	 |  x6+x4+x3+x+1
a33	= a5+a4+a3+a2+a+1	= 111111		= 63	= a96	= a159	= a222	...	 |  x6+x5+x4+x2+1
a34	= a5+a2+1			= 100101		= 37	= a97	= a160	= a223	...	 |  x6+x+1
a35	= a4+1				= 010001		= 17	= a98	= a161	= a224	...	 |  x6+x3+1
a36	= a5+a				= 100010		= 34	= a99	= a162	= a225	...	 |  x3+x+1
a37	= a4+a3+a2+a+1		= 011111		= 31	= a100	= a163	= a226	...	 |  x6+x5+x2+x+1
a38	= a5+a4+a3+a2+a		= 111110		= 62	= a101	= a164	= a227	...	 |  x6+x5+x4+x+1
a39	= a5+a2+a+1			= 100111		= 39	= a102	= a165	= a228	...	 |  x6+x4+x2+x+1
a40	= a4+a2+1			= 010101		= 21	= a103	= a166	= a229	...	 |  x6+x+1
a41	= a5+a3+a			= 101010		= 42	= a104	= a167	= a230	...	 |  x6+x5+x4+x+1
a42	= a3+a2+a+1			= 001111		= 15	= a105	= a168	= a231	...	 |  x2+x+1
a43	= a4+a3+a2+a		= 011110		= 30	= a106	= a169	= a232	...	 |  x6+x5+1
a44	= a5+a4+a3+a2		= 111100		= 60	= a107	= a170	= a233	...	 |  x6+x5+x2+x+1
a45	= a5+a+1			= 100011		= 35	= a108	= a171	= a234	...	 |  x3+x2+1
a46	= a4+a3+a2+1		= 011101		= 29	= a109	= a172	= a235	...	 |  x6+x5+1
a47	= a5+a4+a3+a		= 111010		= 58	= a110	= a173	= a236	...	 |  x6+x5+x3+x2+1
a48	= a5+a3+a2+a+1		= 101111		= 47	= a111	= a174	= a237	...	 |  x6+x5+x4+x2+1
a49	= a2+1				= 000101		= 5		= a112	= a175	= a238	...	 |  x6+x3+1
a50	= a3+a				= 001010		= 10	= a113	= a176	= a239	...	 |  x6+x5+x2+x+1
a51	= a4+a2				= 010100		= 20	= a114	= a177	= a240	...	 |  x6+x4+x2+x+1
a52	= a5+a3				= 101000		= 40	= a115	= a178	= a241	...	 |  x6+x5+x4+x+1
a53	= a3+a+1			= 001011		= 11	= a116	= a179	= a242	...	 |  x6+x5+1
a54	= a4+a2+a			= 010110		= 22	= a117	= a180	= a243	...	 |  x3+x2+1
a55	= a5+a3+a2			= 101100		= 44	= a118	= a181	= a244	...	 |  x6+x5+x3+x2+1
a56	= a+1				= 000011		= 3		= a119	= a182	= a245	...	 |  x6+x3+1
a57	= a2+a				= 000110		= 6		= a120	= a183	= a246	...	 |  x6+x4+x2+x+1
a58	= a3+a2				= 001100		= 12	= a121	= a184	= a247	...	 |  x6+x5+1
a59	= a4+a3				= 011000		= 24	= a122	= a185	= a248	...	 |  x6+x5+x3+x2+1
a60	= a5+a4				= 110000		= 48	= a123	= a186	= a249	...	 |  x6+x4+x2+x+1
a61	= a5+a4+a3+a+1		= 111011		= 59	= a124	= a187	= a250	...	 |  x6+x5+x3+x2+1
a62	= a5+a3+a2+1		= 101101		= 45	= a125	= a188	= a251	...	 |  x6+x5+x3+x2+1
*/