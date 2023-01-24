#include "bch4830.hpp"
using namespace std;

bitset <n> BCH_code_short_t3::generate_data() {
	bitset <n> Data;
	for (int i = 0; i < k; i++) {
		Data[i] = (rand() % 2);
	}
	return Data;
}

void BCH_code_short_t3::read_p() {
// Primitive polynomial of degree 6 - 1011011
	p = 0b1011011;
	primitive_polynomial = p.to_ulong();
	cout << "Primitive polynomial:" << endl << "p(x) = ";
	verbose_polynomial(p);
}

void BCH_code_short_t3::generate_gf() {
	index_of[0] = -1;
	int m = MSB(p);
	for (int i = 0; i < GF; i++) {
		if (i < m) {
			alpha_to[i] = 1 << i;
			index_of[alpha_to[i]] = i;
		} else {
			int previous_alpha = alpha_to[i - 1];
			if (previous_alpha >= 32) {
				alpha_to[i] = (previous_alpha << 1) ^ primitive_polynomial;
			} else {
				alpha_to[i] = previous_alpha << 1;
			}
			index_of[alpha_to[i]] = i;
		}
	}
}

void BCH_code_short_t3::gen_poly() {
/*  
	Compute generator polynomial of BCH code of n(2**6)
*/
	vector <vector <int>>  cycle_cosets;
	set <int> unique_numbers;
	// Generate cycle sets modulo 63
	for (int first_element = 1; first_element <= 31; first_element++) {
		vector <int> coset;
		coset.push_back(first_element);
		int j = -1;
		while (true) {
			j++;
			coset.push_back((coset[j] << 1) % GF);
			//check if element is unique
			auto status = unique_numbers.emplace(coset[j]);
			if (!status.second) {
				break;
			}
			// if the first element in coset will be equal to the next element so that we know if 
			// we made a full circle and its time to push to vector
			if (coset[0] == (coset[j+1] << 1) % GF) {
				cycle_cosets.push_back(coset);
				break;
			}
		}
	}
	// Search for roots 1, 2, ..., t*2 in cycle sets
	for (const auto& coset : cycle_cosets) {
		bool root_found = false;
		for (const auto& alpha : coset) {
			for (int root = 1; root < t*2; root++)
				if (alpha == root) {
					root_found = true;
					break;
				}
			if(root_found) { break;}
		}
		if(root_found) {
		//populate zeros with cosets that have roots 1 to d-1
			zeros_deluxe.push_back(coset);
		}
		if (zeros_deluxe.size() == t) {
			break;
		}
	}

	//calculate first and second minimal polynomial
	vector <int> min_polynomials;
	ULL first_factor, second_factor;
	// multiply all elements from one zero coset and then all elements from second zero coset
	for (const auto& zero_coset : zeros_deluxe) {
		ULL product = 0;
		for (uint i=1; i<zero_coset.size(); i++) {
			if (i == 1) {
				first_factor = alpha_to[zero_coset[i-1]] ^ 2; // (ax + x)
			} else {
				first_factor = product;
			}
			second_factor = alpha_to[zero_coset[i]] ^ 2;
			product = multiply_uint_polynomials(first_factor, second_factor);
		}
		product %= GF+1;
		product ^= primitive_polynomial;
		min_polynomials.push_back(product);
	}

	cout << "This is a (" << n << "," << k << "," << t*2+1 << ") binary BCH code" << endl;

	// Compute generator polynomial by multiplying zeros root polynomials
	uint generator_polynomial = min_polynomials[0];
	for (uint i=1; i<min_polynomials.size(); i++) {
		generator_polynomial =  multiply_uint_polynomials(generator_polynomial, min_polynomials[i]);
	}
	// generator_polynomial = 1100100100111 but the program won't work with it lol so we have to manually set another value
	generator_polynomial_bitset = generator_polynomial;
	reverse_bitset(generator_polynomial_bitset, k-1);
	cout << "g(x) is set to " << print_wihtout_zeros(generator_polynomial_bitset, n-k+1) << endl;
}

bitset <n> BCH_code_short_t3::encode_bch(const bitset <n> &Data) {
/*
	codeword is c(X) = Data(X)*X**(n-k)+ rb(X), data shifted by n-k bits and xored with redundant bits
*/
	// to calculate redundant bits, get the remainder of Data shifted by n-k bits divided by generator polynomial
	auto Shifted_Data = Data<<(n-k);
	auto rb = divide_bitset_polynomials(Shifted_Data, generator_polynomial_bitset).first;
	cout << "DATA       " << Data << endl << "SHFTD DATA " << Shifted_Data << endl << "GEN POL    " << generator_polynomial_bitset << endl << "RED BITS   " << rb << endl;
	// systematic encoding: The message as a suffix, first n-k bits are redundancy, last k bits are message
	bitset <n> Codeword = Shifted_Data ^ rb;
	//check if the generated rb are valid
	auto check = divide_bitset_polynomials(Codeword, generator_polynomial_bitset).first;
	if (check != 0) {
		cout<<"redundant bits are not correct: "<<check<<endl;
	}
	cout << "CODEWORD   " << Codeword << endl;
	cout << endl;
	return Codeword;
}

vector <int> BCH_code_short_t3::calculate_syndromes(const bitset <n> &Received_Codeword, bool &syn_error) {
	vector <int> syndromes(2*t);
	for (int i=1; i<=2*t; i++) {
		syndromes[i] = 0;
		for (int j=0; j<n; j++) {
			if (Received_Codeword[n-j-1] != 0) {
				syndromes[i] ^= alpha_to[(i*j) % GF];
			}
		}
		// convert syndrome from polynomial form to index form
		syndromes[i] = index_of[syndromes[i]];
		if (syndromes[i] != -1) {
			syn_error=true;
		}
	}
	cout<<endl<<"syndromes= ( ";
	for (int i=1; i<=2*t; i++) {
		cout << syndromes[i] << " ";
	}
	cout<<")"<<endl;
	return syndromes;
}

bitset <n> BCH_code_short_t3::decode_bch(const bitset <n> &Received_Codeword) {
	/*
	* Simon Rockliff's implementation of Berlekamp's algorithm.
	* Assume we have received bits in recd[i], i=0..(n-1).
	*
	* Compute the 2*t syndromes by substituting alpha^i into rec(X) and
	* evaluating, storing the syndromes in s[i], i=1..2t (leave s[0] zero) .
	* Then we use the Berlekamp algorithm to find the error location polynomial
	* elp[i].
	*
	* If the degree of the elp is >t, then we cannot correct all the errors, and
	* we have detected an uncorrectable error pattern. We output the information
	* bits uncorrected.
	*
	* If the degree of elp is <=t, we substitute alpha^i , i=1..n into the elp
	* to get the roots, hence the inverse roots, the error location numbers.
	* This step is usually called "Chien's search".
	*
	* If the number of errors located is not equal the degree of the elp, then
	* the decoder assumes that there are more than t errors and cannot correct
	* them, only detect them. We output the information bits uncorrected.
	*/

	bool syn_error;
	// first form the syndromes
	auto s = calculate_syndromes(Received_Codeword, syn_error);
	bitset <n> Decoded_Message = Received_Codeword;
	
	if (syn_error) {
		if (s[3] == (s[1] * 3) % GF) {
			cout << "One error at " << s[1];
			Decoded_Message.flip(n-1-s[1]);
		} else {
			/*
			* Compute the error location polynomial via the Berlekamp
			* iterative algorithm. Following the terminology of Lin and
			* Costello's book :   d[u] is the 'mu'th discrepancy, where
			* u='mu'+1 and 'mu' (the Greek letter!) is the step number
			* ranging from -1 to 2*t (see L&C),  l[u] is the degree of
			* the elp at that step, and u_l[u] is the difference between
			* the step number and the degree of the elp. 
			*
			* params:
			* l --> number of found errors, (the current length of the LFSR) linear feedback shift register
			* d --> discrepancy
			* u_lu
			* err_loc_pol_coeffs[x][y]
			*/
			int err_loc_pol_coeffs[2*t+1][2*t+1] = {{0, 0}};
			int d[2*t+1] = {0}, l[2*t+1] = {0}, u_lu[2*t+1] = {0};
			vector <int> error_locations;

			d[1] = s[1];		/* index form */
			err_loc_pol_coeffs[1][0] = 1;		/* polynomial form */
			int u;
			for (u = 1; l[u + 1] <= t; u++) {
				if (d[u] == -1) {
					l[u + 1] = l[u];
					for (int i = 0; i <= l[u]; i++) {
						err_loc_pol_coeffs[u + 1][i] = err_loc_pol_coeffs[u][i];
						err_loc_pol_coeffs[u][i] = index_of[err_loc_pol_coeffs[u][i]];
					}
				} else {
					/*
					* search for words with greatest u_lu[q] for
					* which d[q]!=0 
					*/
					int q = u - 1;
					while ((d[q] == -1) && (q > 0))
						q--;
					/* have found first non-zero d[q]  */
					if (q > 0) {
						int j = q;
						do {
							j--;
							if ((d[j] != -1) && (u_lu[q] < u_lu[j]))
							q = j;
						} while (j > 0);
					}
					/*
					* have now found q such that d[u]!=0 and
					* u_lu[q] is maximum 
					*/
					/* store degree of new err_loc_pol_coeffs polynomial */
					if (l[u] > l[q] + u - q)
						l[u + 1] = l[u];
					else
						l[u + 1] = l[q] + u - q;

					/* form new err_loc_pol_coeffs(x) */
					for (int i = 0; i < 2*t; i++)
						err_loc_pol_coeffs[u + 1][i] = 0;
					for (int i = 0; i <= l[q]; i++)
						if (err_loc_pol_coeffs[q][i] != -1)
							err_loc_pol_coeffs[u + 1][i + u - q] = 
									alpha_to[(d[u] + GF - d[q] + err_loc_pol_coeffs[q][i]) % GF];
					for (int i = 0; i <= l[u]; i++) {
						err_loc_pol_coeffs[u + 1][i] ^= err_loc_pol_coeffs[u][i];
						err_loc_pol_coeffs[u][i] = index_of[err_loc_pol_coeffs[u][i]];
					}
				}
				u_lu[u + 1] = u - l[u + 1];

				/* form (u+1)th discrepancy */
				/* no discrepancy computed on last iteration */
				if (s[u + 1] != -1) {
					d[u + 1] = alpha_to[s[u + 1]];
				} else {
					d[u + 1] = 0;
				}
				for (int i = 1; i <= l[u + 1]; i++) {
					if ((s[u + 1 - i] != -1) && (err_loc_pol_coeffs[u + 1][i] != 0)) {
						d[u + 1] ^= alpha_to[(s[u + 1 - i] + index_of[err_loc_pol_coeffs[u + 1][i]]) % GF];
					}
				}
				/* put d[u+1] into index form */
				d[u + 1] = index_of[d[u + 1]];	
			}
			
			if (l[u] <= t) {/* Can correct errors */
				/* put err_loc_pol_coeffs into index form */
				cout << "Sigma(x) = ";
				for (int i = 0; i <= l[u]; i++) {
					err_loc_pol_coeffs[u][i] = index_of[err_loc_pol_coeffs[u][i]];
					cout << err_loc_pol_coeffs[u][i] << " ";
				}
				// Find roots of the error location polynomial:
				// Chien search
				cout << endl << "Roots: ";
				for (int i = 1; i <= GF; i++) {
					int q = 1;
					for (int j = 1; j <= l[u]; j++) {
						err_loc_pol_coeffs[u][j] = (err_loc_pol_coeffs[u][j] + j) % GF;
						q ^= alpha_to[err_loc_pol_coeffs[u][j]];
					}
					// Store error location number indices
					if (!q) {
						error_locations.push_back(GF - i);
						cout << error_locations.back() << " ";
					}
				}
				if (error_locations.size() == (unsigned)l[u]) {
				/* no. roots = degree of err_loc_pol_coeffs hence <= t errors */
					for (auto const &error_location : error_locations) {
                        if (error_location <0 || error_location>=n) {
                            cout<<endl<<"Incomplete decoding: errors detected"<<endl;
						    return -1;
                        } else {
                            Decoded_Message.flip((n-1-error_location + n) % n);
                        }
					}
				} else {/* err_loc_pol_coeffs has degree >t hence cannot solve */
					cout<<endl<<"Incomplete decoding: errors detected"<<endl;
                    return -1;
				}
			} else {
				cout<<endl<<"Incomplete decoding: errors detected"<<endl;
                return -1;
			}
		}
	} else {
		cout << "No errors found";
	}
	return Decoded_Message;
}

void BCH_code_short_t3::verbose_polynomial(const bitset <n> &polynomial) { //human readable polynomial format
	int power = MSB(polynomial);
	for (int i=power; i>=0; i--) {
		if (polynomial[i]) {
			if (i != power) {
				cout << " + " ;
			}
			if (i!=0) {
				cout << "x^" << i;
			} else {
				cout << "1";
			}
		}
	}
	cout << endl;
}

template<size_t N>
void BCH_code_short_t3::reverse_bitset(bitset <N> &polynomial, int shift) {
    for(size_t i = 0; i < N/2; ++i) {
    	bool temp_bit = polynomial[i];
    	polynomial[i] = polynomial[N-i-1];
    	polynomial[N-i-1] = temp_bit;
    }
	polynomial >>= shift;
}

template <size_t N>
int BCH_code_short_t3::MSB(const bitset <N> &polynomial) { // Most Significant Bit
	return (N + (GF-N) - countl_zero(polynomial.to_ullong()));
}

pair<bitset <n>, bitset <n>> BCH_code_short_t3::divide_bitset_polynomials(const bitset <n> &dividend, const bitset <n> &divisor) {
	// 	101 0101 0000 0000 | 1 1101 0001
	// _________________|________________
	// 			. . .   |   111 0101 <- quotient
	// 		__________|
	// 		1110 0101 <- remainder
	// division is the same as polynomial modulo polynomial
	bitset <n> quotient, remainder = dividend;;
 	while (MSB(remainder) >= MSB(divisor)) {
		int shift =  MSB(remainder) - MSB(divisor);
		remainder ^= divisor  <<  shift;
		quotient.flip(shift); 
	}
	return {remainder, quotient};
}

uint BCH_code_short_t3::multiply_uint_polynomials(uint mulitplicand, uint multiplicator) {
	uint product = 0;
	while (mulitplicand > 0) {
		if (mulitplicand & 1) {
			product ^= multiplicator;
		}
		multiplicator <<= 1;
		mulitplicand >>= 1;
	}
	return product;
}

bitset <n> BCH_code_short_t3::multiply_bitset_polynomials(const bitset <n> &mulitplicand, const bitset <n> &multiplicator) {
	bitset <n> product;
	for (int i = 0; i < n; ++i) {
		if (multiplicator[i]) {
		product ^= mulitplicand << i;
		}
  	}
	return product;
}


bitset <n> BCH_code_short_t3::user_input(const bitset <n> &Codeword) {
	cout << endl << "Enter the number of errors (choose a number between 1 and 10): " << endl;
	uint numerr;
	cin >> numerr;
	while (!(numerr >= 1 && numerr <= 10)) {
		if (!cin) {
			cout  << "Please enter an integer" << endl;
			cin_clean();
		} else {
			cout << "Wrong number of errors" << endl;
			cin_clean();
		}
		cin >> numerr;
	}
	cout << "Enter the position of errors (reading from left to right choose a number between 0 and "<<n-1<<"): " << endl;
	uint error_position;
	bitset <n> Received_Codeword = Codeword;
	for (uint i = 0; i < numerr; i++) {
		cout << "Position " << i+1 << ":" << endl;
		bool cin_check = false;
		cin >> error_position;
		while (!cin_check) {
			if (!cin) { //we have to check it like that because if cin reads a character instead of integer the cin fails and the 
				cout << "Please enter an integer" << endl;  //check_bad_input function will be passed a value of 0 which would make the previous while loop break
				cin_clean();
				cin >> error_position;
			} else if(!(error_position >= 0 && error_position <= n-1)) {
			cout << "Wrong error position" << endl;
			cin_clean();
			cin >> error_position;
			} else {
			cin_check = true;
			}
		}
		Received_Codeword.flip(n - 1 - error_position);
	}
	return Received_Codeword;
}

void BCH_code_short_t3::cin_clean() {
	cin.clear();
	cin.ignore(numeric_limits<streamsize>::max(), '\n');
}

void BCH_code_short_t3::print_codeword_and_received_codeword(const bitset <n> &Codeword, const bitset <n> &Received_Codeword) {
	cout << "c(x) = "<< Codeword << endl;
	cout << "r(x) = "<< Received_Codeword << endl;
	cout <<"       ";
	for (int i=n-1; i >=0; i--) {
		if (Codeword[i] != Received_Codeword[i]) {
			cout << '^';
		} else {
			cout << ' ';
		}
	}
	cout << endl << "Positions of errors in the received codeword (at ^)" << endl;
}

void BCH_code_short_t3::print_message_and_decoded_message(const bitset <n> &Data, const bitset <n> &Decoded_Data) {
	cout << endl << "Results: " << endl;
	cout << "Original Data  = " << print_wihtout_zeros(Data, k) << endl;
	cout << "Recovered Data = " << print_wihtout_zeros(Decoded_Data>>(n-k), k) << endl;
	// decoding errors: we compare only the Data portion 
	bitset <n> test = Data ^ (Decoded_Data >> (n-k));
	cout<<"                 "<<(test.to_string(' ','^')).substr(n-k)<<endl;
	if (test.count()) {
		cout <<"Position of " << test.count() << " message decoding errors (at ^)\n";
	} else {
		cout <<"Succesful decoding\n";
	}
}

string BCH_code_short_t3::print_wihtout_zeros(const bitset <n> &Polynomial, const uint &Not_Zeros) {
	return (Polynomial.to_string().substr(n - Not_Zeros));
}

int main() {
	char run_program = 'y';
	BCH_code_short_t3 BCH_obj;
	while(run_program == 'y') {
		int seed = time(NULL);
		cout << "Seed used: " << seed << endl;
		srand(seed);
		// Randomly generate message Data
		bitset<n> Data = BCH_obj.generate_data();
		// ecnode message into polynomial
		bitset<n> Codeword = BCH_obj.encode_bch(Data);
		// input errors into codeword
		bitset<n> Received_Codeword = BCH_obj.user_input(Codeword);
		// show codeword and received codeword
		BCH_obj.print_codeword_and_received_codeword(Codeword, Received_Codeword);
		// decode received codeword
		bitset<n> Decoded_Data = BCH_obj.decode_bch(Received_Codeword);
		// print out orignial message and decoded message
        if (Decoded_Data != -1) {
		    BCH_obj.print_message_and_decoded_message(Data, Decoded_Data);
        }
		cout << "Run program again? (y/n)" << endl;
		cin >> run_program;
		while(!cin || (run_program != 'y' && run_program != 'n')) {
			cout << "Run program again? (y/n)\r" << endl;
			BCH_obj.cin_clean();
			cin >> run_program;
		}
	}
    return 0;
}