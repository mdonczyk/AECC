#include "bch6351.hpp"

std::bitset <GF> BCH_code::generate_data() {
	std::bitset <GF> Data;
	for (int i = 0; i < k; i++) {// k = 1
		Data[i] = (rand() % 2);
	}
	return Data;
}

void BCH_code::read_p() { //DONE
// Primitive polynomial of degree 6 - 1011011
	p = 0b1011011;
	prim_polynomial = p.to_ulong();
	std::cout << "Primitive polynomial:" << std::endl << "p(x) = ";
	verbose_polynomial(p);
}

void BCH_code::generate_gf() { //DONE
	uint mask = 1, temp_prim_polynomial = prim_polynomial;
	for (uint i = 0; i < m; i++) { //m = 6
		alpha_to[i] = mask;
		index_of[alpha_to[i]] = i;
		if (temp_prim_polynomial & 1) {
			alpha_to[m] ^= mask;
		}
		temp_prim_polynomial >>= 1;
		mask <<= 1;
	}
	index_of[alpha_to[m]] = m;
	for (uint i = m + 1; i < n + 1; i++) {
		if (alpha_to[i - 1] >= 32) {
			alpha_to[i] = (alpha_to[i - 1] << 1) ^ prim_polynomial; //prim_polynomial = x^6 + x^4 + x^3 + x^1 + x^0 = 1011011 = 91
		}
		else {
			alpha_to[i] = alpha_to[i - 1] << 1;
		}
		index_of[alpha_to[i]] = i;
	}
	index_of[0] = -1;
}

void BCH_code::gen_poly() { //FIXME: second polynomial isn't calculated right, should be 117
/*  
	Compute generator polynomial of BCH code of n = 63, redundancy = 12
*/
	std::vector <std::vector <int>>  cycle_coset;
	std::vector <int> temp_coset_index;
	std::set <int> unique_numbers;
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
			//if current element equals first element then we made a full circle and its time to push to std::vector
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
			for (int root = 1; root < d; root++)
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
		if (roots_found == d - 1) {
			break;
		}
	}
	//calculate first and second minimal polynomial
	std::vector <int> min_polynomials;
	int first_factor, second_factor, product = 0;
	for (const auto& zero_coset : zeros_deluxe) {
		for (int i=1; i<zero_coset.size(); i++) {
			if (i == 1) {
				first_factor = alpha_to[zero_coset[i-1]] ^ 2; // (ax + 1)
			} else {
				first_factor = product;
			}
			second_factor = alpha_to[zero_coset[i]] ^ 2;
			multiply_polynomials(first_factor, second_factor, product);
			product %= n;
		}
		product ^= prim_polynomial;
		min_polynomials.push_back(product);
		product = 0;
	}			
	//FIXME: second polynomial isn't calculated right, should be 117
	min_polynomials[1] = 117;
	std::cout << "This is a (" << n << "," << k << "," << d << ") binary BCH code" << std::endl;
	// Compute generator polynomial by multiplying zeros root polynomials
	for (int i=1; i<min_polynomials.size(); i++) {
		multiply_polynomials(min_polynomials[i-1], min_polynomials[i], generator_polynomial);
	}
	std::bitset <13> generator_polynomial_binary = generator_polynomial;
	std::cout << "g(x) = " << generator_polynomial_binary << std::endl;
}

std::bitset <GF> BCH_code::encode_bch(const std::bitset <GF> &Data) { //DONE multiply message polynomial by generator polynomial
/*
	codeword is c(X) = Data(X)*X**(n-k)+ rb(X)1 n = 63, k = 51;
*/
	std::bitset <GF> Codeword;
	std::bitset <GF> generator_polynomial_binary = generator_polynomial;
	std::bitset <GF> rb = divide_polynomials(Data, generator_polynomial_binary);
	std::cout << "DATA      " << Data << std::endl << "GEN POL   " << generator_polynomial_binary << std::endl << "RED BITS  " << rb << std::endl;
	//shift rb by k and add redundant bits as prefix
	rb <<= k;
	Codeword = Data ^ rb;
	std::cout << "CODEWORD: " << Codeword << std::endl; //Systematic encoding: The message as a suffix
	// first (n-k) bits are Data
	// last k bits are redundancy
	std::cout << "stored c: ";
	for (int i = 0; i < n; i++){
		c[i] = Codeword[n-i-1]; //store codeword
		std::cout << c[i];
	}
	std::cout << std::endl;
	return Codeword;
}

std::bitset <GF> BCH_code::decode_bch(const std::bitset <GF> &Received_Codeword) { //FIXME: rewrite the whole function
/*
	We do not need the Berlekamp algorithm to decode.
	We solve before hand two equations in two variables.
*/
	std::bitset <GF> Decoded_Message = Received_Codeword;
	int i, j, q;
	int elp[3] = {0}, s[5] = {0}, s3;
	int count = 0, syn_error = 0;
	int loc[3] = {0}, err[3] = {0}, reg[3] = {0};
	int	aux;
	// first form the syndromes
	std::cout << std::endl << "s[] = (";
	for (i = 1; i <= 4; i++) {
		s[i] = 0;
		for (j = 0; j < n; j++)
			if (Received_Codeword[j] != 0)
				s[i] ^= alpha_to[(i * j) % n];
		if (s[i] != 0)
			syn_error = 1;	/* set flag if non-zero syndrome 
								NOTE: If only error detection is needed,
								then exit the program here...
							*/
		// convert syndrome from polynomial form to index form  
		s[i] = index_of[s[i]];
		if (i<4)
			std::cout << s[i] << "  ";
		else
			std::cout << s[i] << ")" << std::endl;
	}
	if (syn_error) {	// If there are errors, try to correct them 
		if (s[1] != -1) {
			s3 = (s[1] * 3) % n;
			if ( s[3] == s3 )  // Was it a single error ? 
				{
				std::cout << "One error at " << s[1];
				Decoded_Message[s[1]] = Decoded_Message[s[1]] ^ 1;		// Yes: Correct it 
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
				std::cout << "Sigma(x) = ";
				for (i = 0; i <= 2; i++)
					std::cout << elp[i] << " ";
				std::cout << std::endl << "Roots: ";
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
						std::cout << (i%n) << " ";
					}
				}
				if (count == 2)	
				// no. roots = degree of elp hence 2 errors 
				for (i = 0; i < 2; i++)
						Decoded_Message[loc[i]] = Decoded_Message[loc[i]] ^ 1;
				else	// Cannot solve: Error detection 
					std::cout << std::endl << "Incomplete decoding";
				}
			}
		else if (s[2] != -1) // Error detection 
			std::cout << std::endl << "Incomplete decoding";
	}
	return Decoded_Message;
}

void BCH_code::verbose_polynomial(const auto &polynomial) { //human readable polynomial format //DONE
	uint power = MSB(polynomial);
	unsigned long long mask = pow(2, power);
	unsigned long long init_mask = mask;
	while (mask > 0) {
		if(polynomial.to_ulong() & mask) {
			if (mask == init_mask) {
			std::cout << "x^" << power;
			} else if (mask == 1) {
			std::cout << " + 1";
			} else {
			std::cout << " + x^" << power;
			}
		}
		mask >>= 1;
		power --;
	}
	std::cout << std::endl;
}

uint BCH_code::MSB(const auto &polynomial) { //DONE Most Significant Bit
	return (63 - std::countl_zero(polynomial.to_ullong()));
}

std::bitset <GF> BCH_code::divide_polynomials(const std::bitset <GF> &polynomial1, const std::bitset <GF> &polynomial2) { //DONE
	// 	101 0101 0000 0000 | 1 1101 0001
	// _________________|________________
	// 			. . .   |   111 0101 <- quotient
	// 		__________|
	// 		1110 0101 <- remainder
	// :param polynomial1: 1st polynomial.
	// :param polynomial2: 2nd polynomial.
	// :returns: the quotient and the remainder.
	//std::bitset <GF> quotient, shifter = 1;
	std::bitset <GF> remainder;
	remainder = polynomial1;
	while (MSB(remainder) >= MSB(polynomial2)) {
		int shift =  MSB(remainder) - MSB(polynomial2);
		remainder ^= polynomial2  <<  shift;
		//quotient ^= shifter  <<  shift;
	}
	return remainder;
}

void BCH_code::multiply_polynomials(auto polynomial1, auto polynomial2, auto &product) { //DONE
	while (polynomial1 > 0) {
		if (polynomial1 & 1) {
			product ^= polynomial2;
		}
		polynomial2 <<= 1;
		polynomial1 >>= 1;
	}
}

std::bitset <GF> BCH_code::user_input(const std::bitset <GF> &Codeword) { //DONE
	std::cout << std::endl << "Enter the number of errors (choose a number between 1 and 10): " << std::endl;
	uint numerr;
	std::cin >> numerr;
	while (!(numerr >= 1 && numerr <= 10)) {
		if (!std::cin) {
			std::cout  << "Please enter an integer" << std::endl;
			cin_clean();
		} else {
			std::cout << "Wrong number of errors" << std::endl;
			cin_clean();
		}
		std::cin >> numerr;
	}
	std::cout << "Enter the position of errors (choose a number between 0 and 62): " << std::endl;
	uint error_position;
	unsigned long long shifter = 1;
	std::bitset <GF> Received_Codeword = Codeword;
	for (int i = 0; i < numerr; i++) {
		std::cout << "Position " << i+1 << ":" << std::endl;
		bool cin_check = false;
		std::cin >> error_position;
		while (!cin_check) {
			if (!std::cin) { //we have to check it like that because if std::cin reads a character instead of integer the std::cin fails and the 
				std::cout << "Please enter integer" << std::endl;  //check_bad_input function will be passed a value of 0 which would make the previous while loop break
				cin_clean();
				std::cin >> error_position;
			} else if(!(error_position >= 0 && error_position <= 62)) {
			std::cout << "Wrong error position" << std::endl;
			cin_clean();
			std::cin >> error_position;
			} else {
			cin_check = true;
			}
		}
		errpos.push_back(error_position);
		std::cout << "Received_Codeword before " << Received_Codeword << std::endl;
		Received_Codeword ^= (shifter << (GF-1 - error_position));
		std::cout << "Received_Codeword after  " << Received_Codeword << std::endl;
	}
	return Received_Codeword;
}

void BCH_code::cin_clean() { //DONE
	std::cin.clear();
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

void BCH_code::print_codeword_and_received_codeword(const std::bitset <GF> &Codeword, const std::bitset <GF> &Received_Codeword) {
	bool error_pos_show [63] = {false};
	std::cout << std::endl << "c(x) = "<< Codeword << std::endl;
	std::cout << std::endl << "r(x) = "<< Received_Codeword << std::endl;
	for (int i = 0; i < errpos.size(); i++) { //whats happening here?
		error_pos_show[errpos[i]] = !(error_pos_show[errpos[i]]);
	}
	std::cout << "errpos:";
	for (int i = 0; i < 63; i++) {
		if (error_pos_show[i]) {
			std::cout << '^';
		} else {
			std::cout << ' ';
		}
	}
	std::cout << std::endl;
	errpos.clear();
}

void BCH_code::print_message_and_decoded_message(const std::bitset <GF> &Data, const std::bitset <GF> &Decoded_Data) {
	std::cout << std::endl << "Results: " << std::endl;
	std::cout << "Original Data  = " << Data << std::endl;

	std::cout << std::endl << "Recovered Data = " << Decoded_Data << std::endl;
	// decoding errors: we compare only the Data portion 
	std::cout << std::endl << "                 ";
	for (int i=Data.size(); i>=0; i--) {
		if (Data[i] != Decoded_Data[i]) {
			decerror++;
			std::cout << '^';
		} else {
			std::cout << ' ';
	}
	// for (int i = n - k; i < n; i++) {
	// 	if (Data[n - i-1] != Received_Codeword[i]) {
	// 		decerror++;
	// 		std::cout << '^';
	// 	} else {
	// 		std::cout << ' ';
	// 	}
	}
	if (decerror) {
		std::cout << std::endl << decerror << " Message decoding errors (at ^)\n\n\n";
		decerror = 0;
	} else {
		std::cout << std::endl << "Succesful decoding\n\n\n";
	}
}

int main() { //DONE
	char run_program = 'y';
	BCH_code BCH_obj;
	while(run_program == 'y') {
		int seed = 1669581010;  //time(NULL);
		std::cout << "Seed used: " << seed << std::endl;
		srand(seed);
		// Randomly generate message Data
		auto Data = BCH_obj.generate_data();
		// ecnode message into polynomial
		auto Codeword = BCH_obj.encode_bch(Data);
		// input errors into codeword
		auto Received_Codeword = BCH_obj.user_input(Codeword);
		// show codeword and received codeword
		BCH_obj.print_codeword_and_received_codeword(Codeword, Received_Codeword);
		// decode received codeword
		//auto Decoded_Data = BCH_obj.decode_bch(Received_Codeword);
		// print out orignial message and decoded message
		//BCH_obj.print_message_and_decoded_message(Data, Decoded_Data);
		std::cout << "Run program again? (y/n)" << std::endl;
		std::cin >> run_program;
		while(!std::cin || run_program != 'y' && run_program != 'n') {
			std::cout << "Run program again? (y/n)\r" << std::endl;
			BCH_obj.cin_clean();
			std::cin >> run_program;
		}
	}
    return 0;
}