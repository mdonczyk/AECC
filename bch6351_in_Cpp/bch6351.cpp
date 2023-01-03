#include "bch6351.hpp"
using namespace std;

bitset <GF> BCH_code::generate_data() {
	bitset <GF> Data;
	for (int i = 0; i < k; i++) {// k = 1
		Data[i] = (rand() % 2);
	}
	return Data;
}

void BCH_code::read_p() { //DONE
// Primitive polynomial of degree 6 - 1011011
	p = 0b1011011;
	primitive_polynomial = p.to_ulong();
	cout << "Primitive polynomial:" << endl << "p(x) = ";
	verbose_polynomial(p);
}

void BCH_code::generate_gf() { //DONE
	int mask = 1, temp_primitive_polynomial = primitive_polynomial;
	int m = MSB(p);
	for (int i = 0; i < m; i++) { //m = 6
		alpha_poly_from_index[i] = mask;
		index_of_alpha_from_poly[alpha_poly_from_index[i]] = i;
		if (temp_primitive_polynomial & 1) {
			alpha_poly_from_index[m] ^= mask;
		}
		temp_primitive_polynomial >>= 1;
		mask <<= 1;
	}
	index_of_alpha_from_poly[alpha_poly_from_index[m]] = m;
	for (int i = m + 1; i < n; i++) {
		if (alpha_poly_from_index[i - 1] >= 32) {
			alpha_poly_from_index[i] = (alpha_poly_from_index[i - 1] << 1) ^ primitive_polynomial; //primitive_polynomial = x^6 + x^4 + x^3 + x^1 + x^0 = 1011011 = 91
		}
		else {
			alpha_poly_from_index[i] = alpha_poly_from_index[i - 1] << 1;
		}
		index_of_alpha_from_poly[alpha_poly_from_index[i]] = i;
	}
	index_of_alpha_from_poly[0] = -1;
}

void BCH_code::gen_poly() { //FIXME: second polynomial isn't calculated right, should be 117
/*  
	Compute generator polynomial of BCH code of n = 63, redundancy = 12
*/
	vector <vector <int>>  cycle_coset;
	vector <int> temp_coset_index;
	set <int> unique_numbers;
	// Generate cycle sets modulo 63
	for (int i = 0; i <= 31; i++) {
		int j = -1;
		while (true) {
			j++;
			if (j == 0) {
				temp_coset_index.push_back(i);
			} else {
				temp_coset_index.push_back((temp_coset_index[j - 1]  <<  1 ) % n);
			}
			//check if element is unique
			auto status = unique_numbers.emplace(temp_coset_index[j]);
			if (!status.second) {
				temp_coset_index.clear();
				break;
			}
			//if current element equals first element then we made a full circle and its time to push to vector
			if (temp_coset_index[0] == (temp_coset_index[j]  <<  1) % n) {
				cycle_coset.push_back(temp_coset_index);
				temp_coset_index.clear();
				break;
			}
		}
	}
	int size = 0, roots_found = 0;
	// Search for roots 1, 2, ..., d-1 in cycle sets (d = 5)
	for (const auto& index : cycle_coset) {
		for (const auto& second_index : index) {
			for (int root = 1; root < t*2+1; root++)
				if (root == second_index) {
					size = index.size();
					roots_found++;
				}
		}
		if(size) {
		//populate zeros with cosets that have roots 1 to d-1
			zeros_deluxe.push_back(index);
			size = 0;
		}
		if (roots_found == t*2) {
			break;
		}
	}
	//calculate first and second minimal polynomial
	vector <int> min_polynomials;
	unsigned long long first_factor, second_factor, product = 0;
	for (const auto& zero_coset : zeros_deluxe) {
		for (uint i=1; i<zero_coset.size(); i++) {
			if (i == 1) {
				first_factor = alpha_poly_from_index[zero_coset[i-1]] ^ 2; // (ax + 1)
			} else {
				first_factor = product;
			}
			second_factor = alpha_poly_from_index[zero_coset[i]] ^ 2;
			product = multiply_uint_polynomials(first_factor, second_factor);
		}
		product %= n+1;
		product ^= primitive_polynomial;
		min_polynomials.push_back(product);
		product = 0;
	}
	cout << "This is a (" << n << "," << k << "," << t*2+1 << ") binary BCH code" << endl;
	// Compute generator polynomial by multiplying zeros root polynomials
	uint generator_polynomial = multiply_uint_polynomials(min_polynomials[0], min_polynomials[1]);
	generator_polynomial_bitset = generator_polynomial;
	cout << "g(x) = " << print_wihtout_zeros(generator_polynomial_bitset, n-k+1) << endl;
}

bitset <GF> BCH_code::encode_bch(const bitset <GF> &Data) { //
/*
	codeword is c(X) = Data(X)*X**(n-k)+ rb(X)1, n = 63, k = 51
*/
	// to calculate redundant bits, get the remainder of Data shifted by n-k bits divided by generator polynomial
	auto rb = divide_bitset_polynomials(Data<<(n-k), generator_polynomial_bitset).first;
	cout << "DATA      " << Data << endl << "GEN POL   " << generator_polynomial_bitset << endl << "RED BITS  " << rb << endl;
	// systematic encoding: The message as a suffix, first n-k bits are redundancy, last k bits are message
	bitset <GF> Codeword = Data ^ rb <<k;
	//check if the generated rb are valid
	
	auto check = divide_bitset_polynomials(0b000101111011001110000100100010110111001001111110010010010010100, generator_polynomial_bitset).first;
	if (check != 0) {
		cout<<"rb are not correct: "<<check<<endl;

	}
	Codeword = 0b000101111011001110000100100010110111001001111110010010010010100;
	cout << "CODEWORD: " << Codeword << endl; 
	cout << "stored c: ";
	for (int i = 0; i < n; i++){
		c[i] = Codeword[n-i-1]; // store codeword
		cout << c[i];
	}
	cout << endl;
	return Codeword;
}

vector <int> BCH_code::calculate_syndromes(const bitset <GF> &Received_Codeword, bool &syn_error) {
	vector <int> syndromes(2*t);
	cout<<endl;
	cout << Received_Codeword ;
	for (int i=1; i<=2*t; i++) {
		syndromes[i] = 0;
		cout<<endl;
		for (int j=0; j<n; j++) {
			if (Received_Codeword[n-j-1] != 0) {
				syndromes[i] ^= alpha_poly_from_index[(i*j) % n];
				// cout<<"a["<<(i * j) % n<<"]="<<alpha_poly_from_index[(i * j) % n]<<" ";
				cout<<syndromes[i]<<" ";
			}
		}
		// convert syndrome from polynomial form to index form  

		syndromes[i] = index_of_alpha_from_poly[syndromes[i]];
		if (syndromes[i] != -1) {
			syn_error=true;
		}
	}
	// cout<<endl;
	// for (auto const bit : index_of_alpha_from_poly) {
	// 	cout<<bit<<" ";
	// }
	cout<<endl<<"syndromes= ( ";
	for (int i=1; i<=4; i++) {
		cout << syndromes[i] << " ";
	}
	cout<<")"<<endl;
	return syndromes;
}

bitset <GF> BCH_code::decode_bch(const bitset <GF> &Received_Codeword) { //FIXME: rewrite the whole function
/*
	We do not need the Berlekamp algorithm to decode.
	We solve before hand two equations in two variables.
*/
// syndrome calculation --> error-location block
	bitset <63> Decoded_Message = Received_Codeword;
	int i, j, q;
	int elp[3] = {0}, s[5] = {0}, s3;
	int count = 0;
	bool syn_error;
	int loc[3] = {0}, err[3] = {0}, reg[3] = {0};
	// first form the syndromes
	calculate_syndromes(Received_Codeword, syn_error);
	if (syn_error) {	// If there are errors, try to correct them 
		if (s[1] != -1) {
			s3 = (s[1] * 3) % n;
			if ( s[3] == s3 )  // Was it a single error ? 
				{
				cout << "One error at " << s[1];
				Decoded_Message[n-1-s[1]] = Decoded_Message[n-1-s[1]] ^ 1;		// Yes: Correct it 
				}
			else {				/* Assume two errors occurred and solve
									for the coefficients of sigma(x), the
									error locator polynomail
								*/
				int	aux;
				if (s[3] != -1)
				aux = alpha_poly_from_index[s3] ^ alpha_poly_from_index[s[3]];
				else
				aux = alpha_poly_from_index[s3];

				elp[0] = 0;
				elp[1] = (s[2] - index_of_alpha_from_poly[aux] + n) % n;
				elp[2] = (s[1] - index_of_alpha_from_poly[aux] + n) % n;
				cout << "Sigma(x) = ";
				for (i = 0; i <= 2; i++)
					cout << elp[i] << " ";
				cout << endl << "Roots: ";
				// find roots of the error location polynomial 
				for (i = 1; i <= 2; i++)
					reg[i] = elp[i];
				count = 0;
				for (i = GF; i >= 1; i--) { // Chien search 
					q = 1;
					for (j = 1; j <= 2; j++)
						if (reg[j] != -1) {
							reg[j] = (reg[j] + j) % n;
							q ^= alpha_poly_from_index[reg[j]];
						}
					if (!q) {	// store error location number indices 
						loc[count] = i % n;
						count++;
						cout << (i%n) << " ";
					}
				}
				if (count == 2)	
				// no. roots = degree of elp hence 2 errors 
				for (i = 0; i < 2; i++)
						Decoded_Message[n-1-loc[i]] = Decoded_Message[n-1-loc[i]] ^ 1;
				else	// Cannot solve: Error detection 
					cout << endl << "Incomplete decoding";
				}
			}
		else if (s[2] != -1) // Error detection 
			cout << endl << "Incomplete decoding";
	}
	return Decoded_Message;
}

void BCH_code::verbose_polynomial(const bitset <GF> &polynomial) { //human readable polynomial format //DONE
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

int BCH_code::MSB(const bitset <GF> &polynomial) { //DONE Most Significant Bit
	return (GF - countl_zero(polynomial.to_ullong()));
}

pair<bitset <GF>, bitset <GF>> BCH_code::divide_bitset_polynomials(const bitset <GF> &dividend, const bitset <GF> &divisor) { //DONE
	// 	101 0101 0000 0000 | 1 1101 0001
	// _________________|________________
	// 			. . .   |   111 0101 <- quotient
	// 		__________|
	// 		1110 0101 <- remainder
	// division is the same as polynomial modulo polynomial
	bitset <GF> quotient, remainder = dividend;;
 	while (MSB(remainder) >= MSB(divisor)) {
		int shift =  MSB(remainder) - MSB(divisor);
		remainder ^= divisor  <<  shift;
		quotient ^= 1ULL  <<  shift;
	}
	return {remainder, quotient};
}

uint BCH_code::multiply_uint_polynomials(uint mulitplicand, uint multiplicator) { //DONE
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

bitset <GF> BCH_code::multiply_bitset_polynomials(const bitset <GF> &mulitplicand, const bitset <GF> &multiplicator) {
	bitset <GF> product;
	for (int i = 0; i < n; ++i) {
		if (multiplicator[i]) {
		product ^= mulitplicand << i;
		}
  	}
	return product;
}


bitset <GF> BCH_code::user_input(const bitset <GF> &Codeword) { //DONE
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
	cout << "Enter the position of errors (choose a number between 0 and 62): " << endl;
	uint error_position;
	bitset <GF> Received_Codeword = Codeword;
	for (uint i = 0; i < numerr; i++) {
		cout << "Position " << i+1 << ":" << endl;
		bool cin_check = false;
		cin >> error_position;
		while (!cin_check) {
			if (!cin) { //we have to check it like that because if cin reads a character instead of integer the cin fails and the 
				cout << "Please enter integer" << endl;  //check_bad_input function will be passed a value of 0 which would make the previous while loop break
				cin_clean();
				cin >> error_position;
			} else if(!(error_position >= 0 && error_position <= 62)) {
			cout << "Wrong error position" << endl;
			cin_clean();
			cin >> error_position;
			} else {
			cin_check = true;
			}
		}
		Received_Codeword ^= (1ULL << (n-1 - error_position));
	}
	return Received_Codeword;
}

void BCH_code::cin_clean() { //DONE
	cin.clear();
	cin.ignore(numeric_limits<streamsize>::max(), '\n');
}

void BCH_code::print_codeword_and_received_codeword(const bitset <GF> &Codeword, const bitset <GF> &Received_Codeword) {
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

void BCH_code::print_message_and_decoded_message(const bitset <GF> &Data, const bitset <GF> &Decoded_Data) {
	cout << endl << "Results: " << endl;
	cout << "Original Data  = " << print_wihtout_zeros(Data, k) << endl;
	cout << "Recovered Data = " << print_wihtout_zeros(Decoded_Data, k) << endl;
	// decoding errors: we compare only the Data portion 
	cout << "                 ";
	for (int i=k-1; i >=0; i--) { //compare only message bits
		if (Data[i] != Decoded_Data[i]) {
			decerror++;
			cout << '^';
		} else {
			cout << ' ';
		}
	}
	if (decerror) {
		cout << endl << "Position of " << decerror << " message decoding errors (at ^)\n\n\n";
		decerror = 0;
	} else {
		cout << endl << "Succesful decoding\n\n\n";
	}
}

string BCH_code::print_wihtout_zeros(const bitset <GF> &Polynomial, const uint &Not_Zeros) {
	return (Polynomial.to_string().substr(GF - Not_Zeros));
}

int main() { //DONE
	char run_program = 'y';
	BCH_code BCH_obj;
	while(run_program == 'y') {
		int seed = 1669581011;  //time(NULL);
		cout << "Seed used: " << seed << endl;
		srand(seed);
		// Randomly generate message Data
		bitset<GF> Data = BCH_obj.generate_data();
		// ecnode message into polynomial
		bitset<GF> Codeword = BCH_obj.encode_bch(Data);
		// input errors into codeword
		bitset<GF> Received_Codeword = BCH_obj.user_input(Codeword);
		// show codeword and received codeword
		BCH_obj.print_codeword_and_received_codeword(Codeword, Received_Codeword);
		// decode received codeword
		bitset<GF> Decoded_Data = BCH_obj.decode_bch(Received_Codeword);
		// print out orignial message and decoded message
		BCH_obj.print_message_and_decoded_message(Data, Decoded_Data);
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