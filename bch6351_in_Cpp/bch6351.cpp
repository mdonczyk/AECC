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

int m = 6, n = 63, k = 51, t = 2, d = 5, length = 63, prim_polynomial = 0;
unsigned long long generator_polynomial = 0;
bitset<7> p;
int alpha_to[64] = {0}, index_of[64] = {0}; 
int c[63], recD[63];
bitset <63> Data, recD_deluxe;
vector <int> errpos, zeros, g;
vector <vector <int>> zeros_deluxe;
int numerr, decerror = 0;

class BCH_code{

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

		void encode_bch() { //multiply message polynomial by generator polynomial 
		/* 
		  codeword is c(X) = Data(X)*X**(n-k)+ rb(X)1 length = 63, k = 51;
		*/	
			
			bitset <63> generator_polynomial_binary = generator_polynomial;
			bitset <63> rb (0);
			divide_polynomials(Data, generator_polynomial_binary, rb);

			cout<<"DATA      "<<Data<<endl<<"GEN POL   "<<generator_polynomial_binary<<endl<<"RED BITS  "<<rb<<endl;

			//shift Data by n-k and add redundant bits as suffix
			//rb <<= 51;
			Data <<= (n-k);
			recD_deluxe = Data ^ rb;
			cout<<"CODEWORD: "<<recD_deluxe<<endl; //Systematic encoding: The message as a suffix
			// first (length-k) bits are redundancy 
			// last k bits are Data 
			cout<<"stored c: ";
			for (int i = 0; i < length; i++){
				c[i] = recD_deluxe[length-i-1]; //store codeword
				recD[i] = c[i];
				cout<<c[i];
			}

			cout<<endl;
		}

		void decode_bch() {
		/*
			We do not need the Berlekamp algorithm to decode.
			We solve before hand two equations in two variables.
		*/
			int i, j, q;
			int elp[3] = {0}, s[5] = {0}, s3;
			int count = 0, syn_error = 0;
			int loc[3] = {0}, err[3] = {0}, reg[3] = {0};
			int	aux;
			// first form the syndromes 
			cout<<endl<<"s[] = (";
			for (i = 1; i <= 4; i++) {
				s[i] = 0;
				for (j = 0; j < length; j++)
					if (recD[j] != 0)
						s[i] ^= alpha_to[(i * j) % n];
				if (s[i] != 0)
					syn_error = 1;	/* set flag if non-zero syndrome 
										NOTE: If only error detection is needed,
										then exit the program here...
									*/
				// convert syndrome from polynomial form to index form  
				s[i] = index_of[s[i]];
				if (i<4)
					cout<<s[i]<<"  ";
				else
					cout<<s[i]<<")"<<endl;
			}
			if (syn_error) {	// If there are errors, try to correct them 
				if (s[1] != -1) {
					s3 = (s[1] * 3) % n;
					if ( s[3] == s3 )  // Was it a single error ? 
						{
						cout<<"One error at "<<s[1];
						recD[s[1]] = recD[s[1]] ^ 1;		// Yes: Correct it 
						}
					else {				/* Assume two errors occurred and solve
											for the coefficients of sigma(x), the
											error locator polynomail
										*/
						if (s[3] != -1)
						aux = alpha_to[s3] ^ alpha_to[s[3]];
						else
						aux = alpha_to[s3];

						elp[0] = 0;
						elp[1] = (s[2] - index_of[aux] + n) % n;
						elp[2] = (s[1] - index_of[aux] + n) % n;
						cout<<"Sigma(x) = ";
						for (i = 0; i <= 2; i++)
							cout<<elp[i]<<" ";
						cout<<endl<<"Roots: ";
						// find roots of the error location polynomial 
						for (i = 1; i <= 2; i++)
							reg[i] = elp[i];
						count = 0;
						for (i = 1; i <= 63; i++) { // Chien search 
							q = 1;
							for (j = 1; j <= 2; j++)
								if (reg[j] != -1) {
									reg[j] = (reg[j] + j) % n;
									q ^= alpha_to[reg[j]];
								}
							if (!q) {	// store error location number indices 
								loc[count] = i % n;
								count++;
								cout<<(i%n)<<" ";
							}
						}
						if (count == 2)	
						// no. roots = degree of elp hence 2 errors 
						for (i = 0; i < 2; i++)
								recD[loc[i]] = recD[loc[i]] ^ 1;
						else	// Cannot solve: Error detection 
							cout<<endl<<"Incomplete decoding";
						}
					}
				else if (s[2] != -1) // Error detection 
					cout<<endl<<"Incomplete decoding";
			}
		}


		void user_input(){

			cout<<endl<<"Enter the number of errors (choose a number between 1 and 10): "<<endl;
			cin>>numerr;
			while (!(numerr >= 1 && numerr <= 10)){
				if (!cin){
					cout <<"Please enter an integer"<<endl;
					cin_clean();
				}
				else{
					cout<<"Wrong number of errors"<<endl;
					cin_clean();
				}
				cin>>numerr;
			}
			cout<<"Enter the position of errors (choose a number between 0 and 62): "<<endl;
			int error_position;
			for (int i = 0; i < numerr; i++) {
				cout<<"Position "<<i+1<<":"<<endl;
				bool cin_check = false;
				cin>>error_position;
				while (!cin_check){
					if (!cin){ //we have to check it like that because if cin reads a character instead of integer the cin fails and the 
						cout<<"Please enter integer"<<endl;  //check_bad_input function will be passed a value of 0 which would make the previous while loop break
						cin_clean();
						cin>>error_position;
					} else if(!(error_position >= 0 && error_position <= 62)){ 
					cout<<"Wrong error position"<<endl;
					cin_clean();
					cin>>error_position;
					} else
					cin_check = true;
				}
				errpos.push_back(error_position);
				recD[error_position] = recD[error_position] ^ 1;
			}
		}

		void cin_clean(){
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
		}

	private:
		void read_p(){
		// Primitive polynomial of degree 6 - 1011011
			p = 0b1011011;
			prim_polynomial = p.to_ulong();
			cout<<"Primitive polynomial:"<<endl<<"p(x) = ";
			verbose_polynomial(p);
		}			

		void generate_gf() {
		/*
		block length n = 63 --> 64-1  Gf(2**6)
		generate GF(2**m) from the irreducible polynomial p(X) in p[0]..p[m]
		lookup tables:  index->polynomial form   power_of_alpha[] contains j=alpha**i;
		polynomial form -> index form  alpha[j=alpha**i] = i alpha=2 is the
		primitive element of GF(2**m) 
			
							GF(64) : P(x) = x6 + x4 + x3 + x + 1 = 1011011 = 91

		alpha_to[4] = 16)	         	index_of[4] = 2) 	          	  
-----------------------------------------------------------------------------------------------------------------
index of alpha:							power of alpha						 |  MINIMAL POLYNOMIAL (addition of all elements in a given cyclotomic set)
a0	= 1					= 000001		= 1		= a63	= a126	= a189	...	 |  x+1 bo x - a0 = x xor 1 = x + 1 = 0000011
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

			int mask = 1, temp_prim_polynomial = prim_polynomial;
			for (int i = 0; i < m; i++) { //m = 6
				alpha_to[i] = mask;
				index_of[alpha_to[i]] = i;
				if (temp_prim_polynomial & 1) {
					alpha_to[m] ^= mask;
				}
				temp_prim_polynomial >>=1;
				mask <<= 1;
			}
			index_of[alpha_to[m]] = m;
			for (int i = m + 1; i < n + 1; i++) {
				if (alpha_to[i - 1] >= 32)
					alpha_to[i] = (alpha_to[i - 1] << 1) ^ prim_polynomial; //prim_polynomial = x^6 + x^4 + x^3 + x^1 + x^0 = 1011011 = 91
				else
					alpha_to[i] = alpha_to[i - 1] << 1;
				index_of[alpha_to[i]] = i;
			}
			index_of[0] = -1;
		}

		void divide_polynomials(const bitset <63> &polynomial1, const bitset <63> &polynomial2, bitset <63> &remainder) {
    
			// 	101 0101 0000 0000 | 1 1101 0001
			// _________________|________________
			// 			. . .   |   111 0101 <- quotient
			// 		__________|
			// 		1110 0101 <- remainder
			// :param polynomial1: 1st polynomial.
			// :param polynomial2: 2nd polynomial.
			// :returns: the quotient and the remainder.

			bitset <63> quotient = 0, shifter = 1;
			int shift = 0;
			remainder = polynomial1;
			while (MSB(remainder) >= MSB(polynomial2)) {
				shift =  MSB(remainder) - MSB(polynomial2);
				remainder ^= polynomial2 << shift;
				quotient ^= shifter << shift;
			}
		}

		int MSB(const auto &polynomial) { //Most Significant Bit
			
			return (63 - countl_zero(polynomial.to_ullong()));
		}

		void verbose_polynomial(const auto &polynomial) { //human readable polynomial format

			int power = MSB(polynomial);
			unsigned long long mask = pow(2, power);
			unsigned long long init_mask = mask;
			while (mask > 0) {
				if(polynomial.to_ulong() & mask) {
					if (mask == init_mask) {
					cout<<"x^"<<power;
					} else if (mask == 1){
					cout<<" + 1";
					} else {
					cout<<" + x^"<<power;
					}
				}
				mask >>= 1;
				power --;
			}
			cout<<endl;
		}

		void multiply_polynomials(auto polynomial1, auto polynomial2, auto &product) {

			while (polynomial1 > 0){
						if (polynomial1 & 1){
							product ^= polynomial2;
						}
						polynomial2 <<= 1;
						polynomial1 >>= 1;
					}
		}

		void gen_poly() {
		/*  
			Compute generator polynomial of BCH code of length = 63, redundancy = 12
		*/
			int allnumbers_size_before_insert = 0;
			vector <vector <int>> cycle_coset;
			vector <int> temp_coset_index;
			set <int> allnumbers;

			// Generate cycle sets modulo 63
			for (int i = 0, j; i <= 31; i++) {
				for (j = 0; ; j++) {
					if (j == 0)
						temp_coset_index.push_back(i);
					else
						temp_coset_index.push_back((temp_coset_index[j - 1] << 1 ) % n);
					
					//check if element is unique
					allnumbers_size_before_insert = allnumbers.size();
					allnumbers.insert(temp_coset_index[j]);
					if (allnumbers_size_before_insert == allnumbers.size()) {
						temp_coset_index.clear();
						break;
					}

					//if current element equals first element then we made a full circle and its time to push to vector
					if (temp_coset_index[0] == (temp_coset_index[j] << 1) % n) {
						cycle_coset.push_back(temp_coset_index);
						temp_coset_index.clear();
						break;
					}
				}
			}
			allnumbers.clear();

			int size = 0, roots_found = 0, alpha_zero = 0;
			// Search for roots 1, 2, ..., d-1 in cycle sets (d = 5)

			for (const auto& index : cycle_coset){
				for (const auto& second_index : index){
					for (int root = 1; root < d; root++)
						if (root == second_index){
							size = index.size();
							alpha_zero = index[0];
							roots_found++;
						}
				}
				if(size != 0) {
				//populate zeros with cosets that have roots 1 to d-1
					zeros_deluxe.push_back(index);
					size = 0;
				}
				if (roots_found == d - 1)
					break;
			}

			//calculate first and second minimal polynomial
			vector <int> min_polynomials;
			int first_factor = 0, second_factor = 0, product = 0;
			// vector <int> zeroooes = {3, 6, 5};
			// for (int i=1; i<zeroooes.size(); i++)  {
			// 	if (i == 1) {
			// 		first_factor = alpha_to[zeroooes[i-1]] ^ 2; // (ax + 1)
			// 	} else {
			// 		first_factor = product;
			// 	}
			// 	second_factor = alpha_to[zeroooes[i]] ^ 2;
			// 	// cout<<"first_factor = "<<first_factor<<endl;
			// 	// cout<<"second_factor = "<<second_factor<<endl;	
			// 	multiply_polynomials(first_factor, second_factor, product);
			// }
			// product %= 6;
			// cout<<"product = "<<product<<endl;
			// product ^= 11;
			// cout<<"product mod 11 = "<<product<<endl;
			

			int dividor = 127;
			for (const auto& zero_coset : zeros_deluxe) {
				for (int i=1; i<zero_coset.size(); i++) {
					if (i == 1) {
						first_factor = alpha_to[zero_coset[i-1]] ^ 2; // (ax + 1)
					} else {
						first_factor = product;
					}
					second_factor = alpha_to[zero_coset[i]] ^ 2;
					// cout<<"first_factor = "<<first_factor<<endl;
					// cout<<"second_factor = "<<second_factor<<endl;	
					multiply_polynomials(first_factor, second_factor, product);
					cout<<first_factor<<" "<<second_factor<<" "<<product<<endl;
					product %= dividor;
				}
				product ^= prim_polynomial;
				if (product == 117)
					cout<<"product % "<<dividor<<" = "<<product<<endl;
				cout<<"product % "<<dividor<<" = "<<product<<endl;
				min_polynomials.push_back(product);
				product = 0;
			}			

			//TODO: second polynomial isn't calculated right, should be 117
			min_polynomials[1] = 117;
			
						
			cout<<"This is a ("<<length<<","<<k<<","<<d<<") binary BCH code"<<endl;

			// Compute generator polynomial by multiplying zeros root polynomials
			for (int i=1; i<min_polynomials.size(); i++) {
				multiply_polynomials(min_polynomials[i-1], min_polynomials[i], generator_polynomial);
			}

			bitset <13> generator_polynomial_binary = generator_polynomial;
			cout<<"g(x) = "<<generator_polynomial_binary<<endl;
		}
};

int main() {
	char run_program = 'y';
	bool error_pos_show [63] = {false};
	BCH_code BCH_magic;
	while(run_program == 'y'){
		int seed = 1669581010;  //time(NULL);
		cout<<"Seed used: "<<seed<<endl;
		srand(seed);
		// Randomly generate Data 
		for (int i = 0; i < k; i++) // k = 1
		 	Data[i] = (rand() % 2);

		// ENCODE 
		BCH_magic.encode_bch();

		// ERRORS
		BCH_magic.user_input(); 
		
		cout<<endl<<"c(x) = ";
		for (int i = 0; i < length; i++) {
			cout<<c[i];
		}
		cout<<endl<<"r(x) = ";
		for (int i = 0; i < length; i++)
			cout<<recD[i];

		for (int i = 0; i < errpos.size(); i++){
			error_pos_show[errpos[i]] = !(error_pos_show[errpos[i]]);
		}
		cout<<endl<<"error: ";
		for (int i = 0; i <= 63; i++){
			if (error_pos_show[i])
				cout<<'^';
			else
				cout<<' ';
		}
		for (int i = 0; i < errpos.size(); i++){ //cleanout the array so that in another run we wont have leftover '^'
			error_pos_show[errpos[i]] = false;
		}
		errpos.clear();

		// DECODE
		BCH_magic.decode_bch();
		// print out original and decoded Data
		cout<<endl<<"Results: "<<endl;
		cout<<"Original Data  = ";
		for (int i = k-1; i >= 0; i--)
			cout<<Data[i];
		cout<<endl<<"Recovered Data = ";
		for (int i = length - k; i < length; i++)
			cout<<recD[i];
		// decoding errors: we compare only the Data portion 
		cout<<endl<<"                 ";
		for (int i = length - k; i < length; i++){
			if (Data[length - i-1] != recD[i]){
				decerror++;
				cout<<'^';
			}
			else
				cout<<' ';
		}
		if (decerror){
			cout<<endl<<decerror<<" Message decoding errors (at ^)\n\n\n";
			decerror = 0;
		}
		else
			cout<<endl<<"Succesful decoding\n\n\n";
		cout<<"Run program again? (y/n)"<<endl;
		cin>>run_program;
		while(!cin || run_program != 'y' && run_program != 'n'){
			cout<<"Run program again? (y/n)\r"<<endl;
			BCH_magic.cin_clean();
			cin>>run_program;
		}
	}
    return 0;
}