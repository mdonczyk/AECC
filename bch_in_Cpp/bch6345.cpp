#include "bch6345.hpp"
using namespace std;

// counters:
int uncaught_errors = 0;
int big_errors = 0;
int total_errors = 0;
bool verbose = false;

// stylistic:
string line(n/2, '*');
string dash_line(n/2, '-');

void BCH::read_p() {
// Primitive polynomial of degree 6 - 1011011
	p = 0b1011011;
	primitive_polynomial = p.to_ulong();
	cout << "Primitive polynomial:" << endl << "p(x) = ";
	verbose_polynomial(p);
}

template <size_t N>
int BCH::MSB(const bitset <N> &Polynomial) { // Most Significant Bit
	return (N + (GF-N) - countl_zero(Polynomial.to_ullong()));
}

void BCH::generate_gf() {
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

void BCH::gen_poly() {
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

	cout << "This is a (" + to_string(n) + "," + to_string(k) + "," + to_string(t*2+1 ) + ") binary BCH code" << endl;

	// Compute generator polynomial by multiplying zeros root polynomials
	uint generator_polynomial = min_polynomials[0];
	for (uint i=1; i<min_polynomials.size(); i++) {
		generator_polynomial =  multiply_uint_polynomials(generator_polynomial, min_polynomials[i]);
	}
	// generator_polynomial = 1100100100111 but the program won't work with it lol so we have to manually set another value
	generator_polynomial_bitset = generator_polynomial;
	reverse_bitset(generator_polynomial_bitset, k-1);
	cout << "g(x) is set to " << print_wihtout_zeros(generator_polynomial_bitset, n-k+1) << endl;
	if (!verbose){
		cout<<"Please be patient, calculating..." << endl;
	}
}

bitset <n> BCH_code_long_t3::encode_bch(const bitset <n> &Data) {
/*
	codeword is c(X) = Data(X)*X**(n-k)+ rb(X), data shifted by n-k bits and xored with redundant bits
*/
	// to calculate redundant bits, get the remainder of Data shifted by n-k bits divided by generator polynomial
	auto Shifted_Data = Data<<(n-k);
	auto rb = divide_bitset_polynomials(Shifted_Data, BCH::generator_polynomial_bitset).first;

	// systematic encoding: The message as a suffix, first n-k bits are redundancy, last k bits are message
	bitset <n> Codeword = Shifted_Data ^ rb;
	//check if the generated rb are valid
	auto check = divide_bitset_polynomials(Codeword, BCH::generator_polynomial_bitset).first;
	if (check != 0) {
		cout<<"redundant bits are not correct: "<<check<<endl;
	}
	if (verbose) {
		cout << endl;
	}
	return Codeword;
}

vector <int> BCH_code_long_t3::calculate_syndromes(const bitset <n> &Received_Codeword, bool &Syn_error) {
	vector <int> syndromes(2*t);
	for (int i=1; i<=2*t; i++) {
		syndromes[i] = 0;
		for (int j=0; j<n; j++) {
			if (Received_Codeword[n-j-1] != 0) {
				syndromes[i] ^= BCH::alpha_to[(i*j) % GF];
			}
		}
		// convert syndrome from polynomial form to index form
		syndromes[i] = BCH::index_of[syndromes[i]];
		if (syndromes[i] != -1) {
			Syn_error=true;
		}
	}
	if (verbose) {
		cout<<endl<<"syndromes= ( ";
		for (int i=1; i<=2*t; i++) {
			cout << syndromes[i] << " ";
		}
		cout<<")"<<endl;
	}
	return syndromes;
}

pair<bitset <k>, bool> BCH_code_long_t3::decode_bch(const bitset <n> &Received_Codeword) {
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

	bool syn_error = false;
	// first form the syndromes
	auto s = calculate_syndromes(Received_Codeword, syn_error);
	bitset <n> Decoded_Codeword = Received_Codeword;
	
	if (syn_error) {
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
		* elp[x][y]
		*/
		vector <int> error_locations;
		int elp[100][100] = {{0, 0}};
		int d[100] = {0}, l[100] = {0}, u_lu[100] = {0};
		for (int i = 1; i < t*2; i++) {
			elp[0][i] = -1;	/* index form */
			elp[1][i] = 0;	/* polynomial form */
		}
		u_lu[0] = -1;
		u_lu[1] = 0;
		d[1] = s[1];		/* index form */
		elp[1][0] = 1;		/* polynomial form */
		int u=0;
		do {
			u++;
			if (d[u] == -1) {
				l[u + 1] = l[u];
				for (int i = 0; i <= l[u]; i++) {
					elp[u + 1][i] = elp[u][i];
					elp[u][i] = BCH::index_of[elp[u][i]];
				}
			} else
				/*
				 * search for words with greatest u_lu[q] for
				 * which d[q]!=0 
				 */
			{
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
				/* store degree of new elp polynomial */
				if (l[u] > l[q] + u - q)
					l[u + 1] = l[u];
				else
					l[u + 1] = l[q] + u - q;
 
				/* form new elp(x) */
				for (int i = 0; i < 2*t; i++)
					elp[u + 1][i] = 0;
				for (int i = 0; i <= l[q]; i++)
					if (elp[q][i] != -1)
						elp[u + 1][i + u - q] = 
                                   BCH::alpha_to[(d[u] + GF - d[q] + elp[q][i]) % GF];
				for (int i = 0; i <= l[u]; i++) {
					elp[u + 1][i] ^= elp[u][i];
					elp[u][i] = BCH::index_of[elp[u][i]];
				}
			}
			u_lu[u + 1] = u - l[u + 1];
 
			/* form (u+1)th discrepancy */
			if (u < 2*t) {	
			/* no discrepancy computed on last iteration */
			  if (s[u + 1] != -1)
			    d[u + 1] = BCH::alpha_to[s[u + 1]];
			  else
			    d[u + 1] = 0;
			    for (int i = 1; i <= l[u + 1]; i++)
			      if ((s[u + 1 - i] != -1) && (elp[u + 1][i] != 0))
			        d[u + 1] ^= BCH::alpha_to[(s[u + 1 - i] 
			                      + BCH::index_of[elp[u + 1][i]]) % GF];
			  /* put d[u+1] into index form */
			  d[u + 1] = BCH::index_of[d[u + 1]];	
			}
		} while ((u < 2*t) && (l[u + 1] <= t));
		
		u++;
		
		// for (int it = 0; it<2*t+1; it++) {
		// 	cout << endl << it << ": ";
		// 	for (int iter = 0; iter<2*t+1; iter++) {
		// 		cout << elp[it][iter] << " ";
		// 	}
		// }

		if (l[u] <= t) {/* Can correct errors */
			/* put elp into index form */
			if (verbose) {
				cout << "Sigma(x) = ";
			}
			for (int i = 0; i <= l[u]; i++) {
				elp[u][i] = BCH::index_of[elp[u][i]];
				if (verbose) {
					cout << elp[u][i] << " ";
				}
			}
			// Find roots of the error location polynomial:
			// Chien search
			if (verbose) {
				cout << endl << "Roots: ";
			}
			for (int i = 1; i <= GF; i++) {
				int q = 1;
				for (int j = 1; j <= l[u]; j++) {
					elp[u][j] = (elp[u][j] + j) % GF;
					q ^= BCH::alpha_to[elp[u][j]];
				}
				// Store error location number indices
				if (!q) {
					error_locations.push_back(n - i + (GF-n));
					if (verbose) {
						cout << error_locations.back() << " ";
					}
				}
			}
			if (error_locations.size() == (unsigned)l[u]) {
			/* no. roots = degree of elp hence <= t errors */
				for (auto const &error_location : error_locations) {
					auto err_loc = (n-1-error_location + n) % n;
					if (err_loc < 0 || err_loc >= n) {
						if (verbose) {
							cout<<endl<<"Incomplete decoding: errors detected"<<endl;
						}
						return {0, false};
					} else {
						Decoded_Codeword.flip(err_loc);
					}
				}
			} else {/* elp has degree >t hence cannot solve */
				if (verbose) {
					cout<<endl<<"Incomplete decoding: errors detected"<<endl;
				}
				return {0, false};
			}
		} else {
			if (verbose) {
				cout<<endl<<"Incomplete decoding: errors detected"<<endl;
			}
			return {0, false};
		}
	} else {
		if (verbose) {
			cout << "No errors found";
		}
	}
	bitset<k> Decoded_Message(Decoded_Codeword.to_string().substr(0, k));
	return {Decoded_Message, true};
}

void BCH::verbose_polynomial(const bitset <n> &Polynomial) { //human readable polynomial format
	int power = MSB(Polynomial);
	for (int i=power; i>=0; i--) {
		if (Polynomial[i]) {
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
void BCH::reverse_bitset(bitset <N> &Polynomial, int Shift) {
    for(size_t i = 0; i < N/2; ++i) {
    	bool temp_bit = Polynomial[i];
    	Polynomial[i] = Polynomial[N-i-1];
    	Polynomial[N-i-1] = temp_bit;
    }
	Polynomial >>= Shift;
}

pair<bitset <n>, bitset <n>> BCH_code_long_t3::divide_bitset_polynomials(const bitset <n> &Dividend, const bitset <n> &Divisor) {
	// 	101 0101 0000 0000 | 1 1101 0001
	// _________________|________________
	// 			. . .   |   111 0101 <- quotient
	// 		__________|
	// 		1110 0101 <- remainder
	// division is the same as polynomial modulo polynomial
	bitset <n> quotient, remainder = Dividend;
 	while (BCH::MSB(remainder) >= BCH::MSB(Divisor)) {
		int shift =  BCH::MSB(remainder) - BCH::MSB(Divisor);
		remainder ^= Divisor  <<  shift;
		quotient.flip(shift); 
	}
	return {remainder, quotient};
}

uint BCH::multiply_uint_polynomials(uint Mulitplicand, uint Multiplicator) {
	uint product = 0;
	while (Mulitplicand > 0) {
		if (Mulitplicand & 1) {
			product ^= Multiplicator;
		}
		Multiplicator <<= 1;
		Mulitplicand >>= 1;
	}
	return product;
}

bitset <n> BCH_code_long_t3::multiply_bitset_polynomials(const bitset <n> &Mulitplicand, const bitset <n> &Multiplicator) {
	bitset <n> product;
	for (int i = 0; i < n; ++i) {
		if (Multiplicator[i]) {
		product ^= Mulitplicand << i;
		}
  	}
	return product;
}

void BCH_code_long_t3::print_codeword_and_received_codeword(const bitset <n> &Codeword, const bitset <n> &Received_Codeword) {
	if (verbose) {
		cout << "c(x) = "<< Codeword << endl;
		cout << "r(x) = "<< Received_Codeword << endl;
	}
	bitset <n> bit_difference = Codeword ^ Received_Codeword;
	if (bit_difference.count()) {
		if (verbose) {
			cout<<"       "<<(bit_difference.to_string(' ','^'))<<endl;
			cout << "Positions of " << bit_difference.count() << " errors in the received codeword at ^.";
		}
		total_errors += bit_difference.count();
		if (bit_difference.count() > t) {
			big_errors++;
		}
		if (verbose) {
			if (bit_difference != 0) {
				cout<< endl << "Positions of errors: ";
			}
			while (bit_difference != 0) {
				int temp = countl_zero(bit_difference.to_ullong());
				cout << temp -(GF - n) -1 << " ";
				bit_difference.flip(GF-temp);
			}
		}
	} else {
		if (verbose) {
			cout << "No errors in the received codeword.";
		}
	}
}

void BCH_code_long_t3::print_message_and_decoded_message(const bitset <n> &Data, const bitset <k> &Decoded_Data) {
	if (verbose) {
		cout << endl << dash_line << "Results" << dash_line << endl;
	}
	bitset<k> Original_Data(Data.to_string().substr(n-k, n));
	if (verbose) {
		cout << "Original Data  = " << Original_Data << endl;
		cout << "Recovered Data = " << Decoded_Data << endl;
	}
	// decoding errors: we compare only the Data portion 
	bitset <k> bit_difference = Original_Data ^ Decoded_Data;
	if (bit_difference.count()) {
		if (verbose) {
			cout << "                 " << (bit_difference.to_string(' ','^')) << endl;
			cout << "Position of " << bit_difference.count() << " message decoding errors at ^." << endl;
		}
		uncaught_errors++;
	} else {
		if (verbose) {
			cout <<"Succesful decoding." << endl;
		}
	}
}

template <size_t N>
string BCH::print_wihtout_zeros(const bitset <N> &Polynomial, const uint &Not_Zeros) {
	return (Polynomial.to_string().substr(N - Not_Zeros));
}

bitset <n> BCH_code_long_t3::introduce_errors(const bitset <n> &Codeword, const int &Probability) {
	bitset <n> modified_codeword = Codeword;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> d(0, Probability);
    for (int i = n - 1; i >= 0; i--) {
        int randnum = d(gen);
        if (randnum == 0) {
            modified_codeword.flip(i);
        }
    }
	return modified_codeword;
}

vector <bitset <k>> bits_to_bitsets (const vector <bool> &Buffer_bits) {
    // put each bit from bool vector to a bitset and then the full bitset to a vector:
    vector <bitset <k>> vector_of_bitsets;
    bitset <k> bits;
    int it = 0;
    for (uint i=0; i<Buffer_bits.size(); i++) {
        bits[it] = Buffer_bits[i];
        it++;
        if (it == k || i == Buffer_bits.size()-1) {
            vector_of_bitsets.push_back(bits);
            it = 0;
            bits = 0;
        }
    }
    return vector_of_bitsets;
}

vector <bool> bitset_to_bits (const vector <bitset <k>> &Vector_of_bitsets) {
// put all bits from bitset into a bool vector:
    vector <bool> recovered_bits;
    for (auto const &bits : Vector_of_bitsets) {
        for (int j = 0; j < k; j++) {
            recovered_bits.push_back(bits[j]);
        }
    }
    return recovered_bits;
}


vector <bool> bytes_to_bits (const vector <unsigned char> &Buffer) {
    // put each bit of the char stream into a bool vector:
    vector <bool> buffer_bits;
     for (int i = 0; i <Buffer.size(); i++) {
        int mask = 0;
            for (int j=0; j<8; j++) {
            bool is_set = Buffer[i] & (1 << mask);
            buffer_bits.push_back(is_set);
            mask ++;
        }
    }
    return buffer_bits;
}

vector <char> bits_to_bytes (const vector <bool> &bits, const int fileSize) {
// put bits from buffer_bits into vector of chars:
    int it = 0;
    vector <char> charstream;
    for (int i = 0; i <fileSize; i++) {
        char temp_char = 0;
            for (int j=0; j<8; j++) {
            temp_char ^= (bits[it] << j);
            it ++;
        }
        charstream.push_back(temp_char);
    }
    return charstream;
}

int main(int argc, char* argv[]) {
	BCH::initialize_BCH();
	ifstream original_image(argv[1], ios::binary);
	int error_probability = stoi(argv[2]);
	if (argc == 4) {
		if (strcasecmp(argv[3], "verbose") == 0) {
			verbose = true;
		}
	}
	ofstream image_with_errors, image_fixed;
    remove("images/image_with_errors_BCH_code_short_t3.bmp");
	remove("images/image_fixed_BCH_code_short_t3.bmp");
    image_with_errors.open ("images/image_with_errors_BCH_code_short_t3.bmp", ios::out | ios::app | ios::binary);
    image_fixed.open ("images/image_fixed_BCH_code_short_t3.bmp", ios::out | ios::app | ios::binary);

    // read the bytes to buffer:
    std::vector<unsigned char> buffer(std::istreambuf_iterator<char>(original_image), {});

	// put each bit of the char stream into a bool vector:
    vector <bool> buffer_bits = bytes_to_bits(buffer);

	// put each bit from bool vector to a bitset and then the full bitset to a vector:
    vector <bitset <k>> vector_of_bitsets = bits_to_bitsets(buffer_bits);

	vector <bitset <k>> vector_of_modified_bitsets(vector_of_bitsets.size());
	vector <bitset <k>> vector_of_recovered_bitsets(vector_of_bitsets.size());

	// create vector of all objects needed for generated bitsets
	vector <BCH_code_long_t3> BCH_objects(vector_of_bitsets.size());

	// more counters:
	int success = 0;
	int failure = 0;

	// start the clock:
	auto start = chrono::high_resolution_clock::now();

	for (int i = 0; i < vector_of_bitsets.size(); i++) {
		if (verbose) {
			cout << line << " START " << line;
		}

		// put generataed bitset as data to object
		BCH_objects[i].Data = bitset<n>(vector_of_bitsets[i].to_string());
	
		// encode message into polynomial
		BCH_objects[i].Codeword = BCH_objects[i].encode_bch(BCH_objects[i].Data);

		// introduce errors into buffer bits
		BCH_objects[i].Received_Codeword = BCH_objects[i].introduce_errors(BCH_objects[i].Codeword, error_probability);

		// show codeword and received codeword
		BCH_objects[i].print_codeword_and_received_codeword(BCH_objects[i].Codeword, BCH_objects[i].Received_Codeword);

		// decode received codeword
		BCH_objects[i].Decoded_Data = BCH_objects[i].decode_bch(BCH_objects[i].Received_Codeword);

		// print out orignial message and decoded message
        if (BCH_objects[i].Decoded_Data.second) {
		    BCH_objects[i].print_message_and_decoded_message(BCH_objects[i].Data, BCH_objects[i].Decoded_Data.first);
			success++;
        } else {
			failure++;
		}
		BCH_objects[i].Received_Data = bitset <k>(BCH_objects[i].Received_Codeword.to_string().substr(0, k));

		// put Received_Data and Decoded_Data to vectors:
		vector_of_modified_bitsets[i] = BCH_objects[i].Received_Data;
		vector_of_recovered_bitsets[i] = BCH_objects[i].Decoded_Data.first;
		
		if (verbose) {
			cout << line << " STOP *" << line << endl << endl;
		}
	}

	// end time clock and show the time:
	auto stop = chrono::high_resolution_clock::now();
	auto dur = stop - start;
	auto duration = std::chrono::duration_cast<std::chrono::duration<float>>(dur);

	// put all bits from bitset into a bool vector:
	vector <bool> modified_bits = bitset_to_bits(vector_of_modified_bitsets);
	vector <bool> recovered_bits = bitset_to_bits(vector_of_recovered_bitsets);

	// put bits from buffer_bits into vector of chars:
	vector <char> modified_charstream = bits_to_bytes(modified_bits, buffer.size());
	vector <char> recovered_charstream = bits_to_bytes(recovered_bits, buffer.size());

    // write all header bytes to new files without modification:
    for (int i = 0; i < HEADER_BYTES; i++) {
        image_with_errors << buffer[i];
		image_fixed << buffer[i];
    }

	// put bits with erros and corrected bits into appropriate files:
	for (int i = HEADER_BYTES; i < modified_charstream.size(); i++) {
		image_with_errors << modified_charstream[i];
		image_fixed << recovered_charstream[i];
	}

	// compare bits from start and those after modification and decoding:
	int difference_count = 0;
	for (int i = 0; i < buffer_bits.size(); i++) {
		if (buffer_bits[i] != recovered_bits[i]) {
			difference_count++;
		}
	}

	// show information about decoding
	cout << endl << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Code used " << setw(20) << "| (" + to_string(n) + "," + to_string(k) + "," + to_string(t*2+1 ) + ")" << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Original image used " << setw(20) << "| " + string(argv[1]) << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Given error probability " << setw(20) << "| 1 in " + to_string(error_probability) << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Real error probability " << setw(20) << "| 1 in " + to_string((vector_of_bitsets.size()*n)/total_errors) << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of all bits " << setw(20) << "| " + to_string(vector_of_bitsets.size() * k) << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of all bits + redundant bits " << setw(20) << "| " + to_string(vector_of_bitsets.size() * n) << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of all generated errors " << setw(20) << "| " + to_string(total_errors) << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of all message polynomials " << setw(20) << "| " + to_string(vector_of_bitsets.size()) << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of successful decodings " << setw(20) << "| " + to_string(success) << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of decoding errors " << setw(20) << "| " + to_string(failure) << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of uncaught decoding errors " << setw(20) << "| " + to_string(uncaught_errors) << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of over t errors in codeword " << setw(20) << "| " + to_string(big_errors) << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Final bit difference " << setw(20) << "| " + to_string(difference_count) << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Encoding and decoding time " << setw(20) << "| " + to_string(duration.count()).substr(0, to_string(duration.count()).length()-3) + " seconds " << "|" << endl;
	cout << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	
	cout << endl << "Results have been written to BCH_logs.txt file" << endl;
	cout << endl << "To view made images on linux use \"feh -F -Z --force-aliasing -d " <<  argv[1] << " images/image_with_errors_BCH_code_short_t3.bmp images/image_fixed_BCH_code_short_t3.bmp\"" << endl << endl;

	// write all the infromation to BCH_logs.txt file
	ofstream BCH_logs;
	BCH_logs.open ("BCH_logs.txt", ios::out | ios::app);
	BCH_logs << endl << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Code used " << setw(20) << "| (" + to_string(n) + "," + to_string(k) + "," + to_string(t*2+1 ) + ")" << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Original image used " << setw(20) << "| " + string(argv[1]) << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Given error probability " << setw(20) << "| 1 in " + to_string(error_probability) << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Real error probability " << setw(20) << "| 1 in " + to_string((vector_of_bitsets.size()*n)/total_errors) << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Number of all bits " << setw(20) << "| " + to_string(vector_of_bitsets.size() * k) << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Number of all bits + redundant bits " << setw(20) << "| " + to_string(vector_of_bitsets.size() * n) << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Number of all generated errors " << setw(20) << "| " + to_string(total_errors) << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Number of all message polynomials " << setw(20) << "| " + to_string(vector_of_bitsets.size()) << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Number of successful decodings " << setw(20) << "| " + to_string(success) << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Number of decoding errors " << setw(20) << "| " + to_string(failure) << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Number of uncaught decoding errors " << setw(20) << "| " + to_string(uncaught_errors) << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Number of over t errors in codeword " << setw(20) << "| " + to_string(big_errors) << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Final bit difference " << setw(20) << "| " + to_string(difference_count) << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	BCH_logs << "|" << setw(37) <<" Encoding and decoding time " << setw(20) << "| " + to_string(duration.count()).substr(0, to_string(duration.count()).length()-3) + " seconds " << "|" << endl;
	BCH_logs << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;

    // close opened files:
    original_image.close();
    image_with_errors.close();
	image_fixed.close();
	BCH_logs.close();
    return 0;
}
