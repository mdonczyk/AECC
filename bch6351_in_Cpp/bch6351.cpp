#include "bch6351.hpp"


void BCH_code::read_p() {
// Primitive polynomial of degree 6 - 1011011
	p = 0b1011011;
	prim_polynomial = p.to_ulong();
	std::cout<<"Primitive polynomial:"<<std::endl<<"p(x) = ";
	verbose_polynomial(p);
}

void BCH_code::generate_gf() {
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

void BCH_code::gen_poly() {
/*  
	Compute generator polynomial of BCH code of n = 63, redundancy = 12
*/
	int allnumbers_size_before_insert = 0;
	std::vector <std::vector <int>> cycle_coset;
	std::vector <int> temp_coset_index;
	std::set <int> allnumbers;
	// Generate cycle sets modulo 63
	for (int i = 0; i <= 31; i++) {
		int j = -1;
		while (true) {
			j++;
			if (j == 0) {
				temp_coset_index.push_back(i);
			} else {
				temp_coset_index.push_back((temp_coset_index[j - 1] << 1 ) % n);
			}
			//check if element is unique
			auto status = allnumbers.emplace(temp_coset_index[j]);
			if (!status.second) {
				temp_coset_index.clear();
				break;
			}
			//if current element equals first element then we made a full circle and its time to push to std::vector
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
	for (const auto& index : cycle_coset) {
		for (const auto& second_index : index) {
			for (int root = 1; root < d; root++)
				if (root == second_index) {
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
		if (roots_found == d - 1) {
			break;
		}
	}
	//calculate first and second minimal polynomial
	std::vector <int> min_polynomials;
	int first_factor = 0, second_factor = 0, product = 0;
	// std::vector <int> zeroooes = {3, 6, 5};
	// for (int i=1; i<zeroooes.size(); i++)  {
	// 	if (i == 1) {
	// 		first_factor = alpha_to[zeroooes[i-1]] ^ 2; // (ax + 1)
	// 	} else {
	// 		first_factor = product;
	// 	}
	// 	second_factor = alpha_to[zeroooes[i]] ^ 2;
	// 	// std::cout<<"first_factor = "<<first_factor<<std::endl;
	// 	// std::cout<<"second_factor = "<<second_factor<<std::endl;	
	// 	multiply_polynomials(first_factor, second_factor, product);
	// }
	// product %= 6;
	// std::cout<<"product = "<<product<<std::endl;
	// product ^= 11;
	// std::cout<<"product mod 11 = "<<product<<std::endl;
	int dividor = n;
	for (const auto& zero_coset : zeros_deluxe) {
		for (int i=1; i<zero_coset.size(); i++) {
			if (i == 1) {
				first_factor = alpha_to[zero_coset[i-1]] ^ 2; // (ax + 1)
			} else {
				first_factor = product;
			}
			second_factor = alpha_to[zero_coset[i]] ^ 2;
			// std::cout<<"first_factor = "<<first_factor<<std::endl;
			// std::cout<<"second_factor = "<<second_factor<<std::endl;	
			multiply_polynomials(first_factor, second_factor, product);
			std::cout<<first_factor<<" "<<second_factor<<" "<<product<<std::endl;
			product %= dividor;
		}
		product ^= prim_polynomial;
		if (product == 117) {
			std::cout<<"product % "<<dividor<<" = "<<product<<std::endl;
		}
		std::cout<<"product % "<<dividor<<" = "<<product<<std::endl;
		min_polynomials.push_back(product);
		product = 0;
	}			
	//TODO: second polynomial isn't calculated right, should be 117
	min_polynomials[1] = 117;
	std::cout<<"This is a ("<<n<<","<<k<<","<<d<<") binary BCH code"<<std::endl;
	// Compute generator polynomial by multiplying zeros root polynomials
	for (int i=1; i<min_polynomials.size(); i++) {
		multiply_polynomials(min_polynomials[i-1], min_polynomials[i], generator_polynomial);
	}
	std::bitset <13> generator_polynomial_binary = generator_polynomial;
	std::cout<<"g(x) = "<<generator_polynomial_binary<<std::endl;
}

void BCH_code::encode_bch() { //multiply message polynomial by generator polynomial 
/*
	codeword is c(X) = Data(X)*X**(n-k)+ rb(X)1 n = 63, k = 51;
*/
	std::bitset <63> generator_polynomial_binary = generator_polynomial;
	std::bitset <63> rb (0);
	divide_polynomials(Data, generator_polynomial_binary, rb);
	std::cout<<"DATA      "<<Data<<std::endl<<"GEN POL   "<<generator_polynomial_binary<<std::endl<<"RED BITS  "<<rb<<std::endl;
	//shift rb by k and add redundant bits as prefix
	rb <<= k;
	recD_deluxe = Data ^ rb;
	std::cout<<"CODEWORD: "<<recD_deluxe<<std::endl; //Systematic encoding: The message as a suffix
	// first (n-k) bits are Data
	// last k bits are redundancy
	std::cout<<"stored c: ";
	for (int i = 0; i < n; i++){
		c[i] = recD_deluxe[n-i-1]; //store codeword
		recD[i] = c[i];
		std::cout<<c[i];
	}
	std::cout<<std::endl;
}

void BCH_code::decode_bch() {
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
	std::cout<<std::endl<<"s[] = (";
	for (i = 1; i <= 4; i++) {
		s[i] = 0;
		for (j = 0; j < n; j++)
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
			std::cout<<s[i]<<"  ";
		else
			std::cout<<s[i]<<")"<<std::endl;
	}
	if (syn_error) {	// If there are errors, try to correct them 
		if (s[1] != -1) {
			s3 = (s[1] * 3) % n;
			if ( s[3] == s3 )  // Was it a single error ? 
				{
				std::cout<<"One error at "<<s[1];
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
				std::cout<<"Sigma(x) = ";
				for (i = 0; i <= 2; i++)
					std::cout<<elp[i]<<" ";
				std::cout<<std::endl<<"Roots: ";
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
						std::cout<<(i%n)<<" ";
					}
				}
				if (count == 2)	
				// no. roots = degree of elp hence 2 errors 
				for (i = 0; i < 2; i++)
						recD[loc[i]] = recD[loc[i]] ^ 1;
				else	// Cannot solve: Error detection 
					std::cout<<std::endl<<"Incomplete decoding";
				}
			}
		else if (s[2] != -1) // Error detection 
			std::cout<<std::endl<<"Incomplete decoding";
	}
}

void BCH_code::verbose_polynomial(const auto &polynomial) { //human readable polynomial format
	int power = MSB(polynomial);
	unsigned long long mask = pow(2, power);
	unsigned long long init_mask = mask;
	while (mask > 0) {
		if(polynomial.to_ulong() & mask) {
			if (mask == init_mask) {
			std::cout<<"x^"<<power;
			} else if (mask == 1) {
			std::cout<<" + 1";
			} else {
			std::cout<<" + x^"<<power;
			}
		}
		mask >>= 1;
		power --;
	}
	std::cout<<std::endl;
}

int BCH_code::MSB(const auto &polynomial) { //Most Significant Bit
	return (63 - std::countl_zero(polynomial.to_ullong()));
}

void BCH_code::divide_polynomials(const std::bitset <63> &polynomial1, const std::bitset <63> &polynomial2, std::bitset <63> &remainder) {
	// 	101 0101 0000 0000 | 1 1101 0001
	// _________________|________________
	// 			. . .   |   111 0101 <- quotient
	// 		__________|
	// 		1110 0101 <- remainder
	// :param polynomial1: 1st polynomial.
	// :param polynomial2: 2nd polynomial.
	// :returns: the quotient and the remainder.
	std::bitset <63> quotient = 0, shifter = 1;
	int shift = 0;
	remainder = polynomial1;
	while (MSB(remainder) >= MSB(polynomial2)) {
		shift =  MSB(remainder) - MSB(polynomial2);
		remainder ^= polynomial2 << shift;
		quotient ^= shifter << shift;
	}
}

void BCH_code::multiply_polynomials(auto polynomial1, auto polynomial2, auto &product) {
	while (polynomial1 > 0) {
				if (polynomial1 & 1) {
					product ^= polynomial2;
				}
				polynomial2 <<= 1;
				polynomial1 >>= 1;
			}
}

void BCH_code::user_input() {
	std::cout<<std::endl<<"Enter the number of errors (choose a number between 1 and 10): "<<std::endl;
	std::cin>>numerr;
	while (!(numerr >= 1 && numerr <= 10)) {
		if (!std::cin) {
			std::cout <<"Please enter an integer"<<std::endl;
			cin_clean();
		} else {
			std::cout<<"Wrong number of errors"<<std::endl;
			cin_clean();
		}
		std::cin>>numerr;
	}
	std::cout<<"Enter the position of errors (choose a number between 0 and 62): "<<std::endl;
	int error_position;
	for (int i = 0; i < numerr; i++) {
		std::cout<<"Position "<<i+1<<":"<<std::endl;
		bool cin_check = false;
		std::cin>>error_position;
		while (!cin_check) {
			if (!std::cin) { //we have to check it like that because if std::cin reads a character instead of integer the std::cin fails and the 
				std::cout<<"Please enter integer"<<std::endl;  //check_bad_input function will be passed a value of 0 which would make the previous while loop break
				cin_clean();
				std::cin>>error_position;
			} else if(!(error_position >= 0 && error_position <= 62)) {
			std::cout<<"Wrong error position"<<std::endl;
			cin_clean();
			std::cin>>error_position;
			} else {
			cin_check = true;
			}
		}
		errpos.push_back(error_position);
		recD[error_position] = recD[error_position] ^ 1;
	}
}

void BCH_code::cin_clean() {
	std::cin.clear();
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

int main() {
	char run_program = 'y';
	bool error_pos_show [63] = {false};
	BCH_code BCH_object;
	while(run_program == 'y') {
		int seed = 1669581010;  //time(NULL);
		std::cout<<"Seed used: "<<seed<<std::endl;
		srand(seed);
		// Randomly generate Data 
		for (int i = 0; i < BCH_object.k; i++) {// k = 1
		 	BCH_object.Data[i] = (rand() % 2);
		}
		// ENCODE 
		BCH_object.encode_bch();
		// ERRORS
		BCH_object.user_input(); 
		std::cout<<std::endl<<"c(x) = ";
		for (int i = 0; i < BCH_object.n; i++) {
			std::cout<<BCH_object.c[i];
		}
		std::cout<<std::endl<<"r(x) = ";
		for (int i = 0; i < BCH_object.n; i++) {
			std::cout<<BCH_object.recD[i];
		}
		for (int i = 0; i < BCH_object.errpos.size(); i++) {
			error_pos_show[BCH_object.errpos[i]] = !(error_pos_show[BCH_object.errpos[i]]);
		}
		std::cout<<std::endl<<"error: ";
		for (int i = 0; i <= 63; i++) {
			if (error_pos_show[i]) {
				std::cout<<'^';
			} else {
				std::cout<<' ';
			}
		}
		for (int i = 0; i < BCH_object.errpos.size(); i++) { //cleanout the array so that in another run we wont have leftover '^'
			error_pos_show[BCH_object.errpos[i]] = false;
		}
		BCH_object.errpos.clear();
		// DECODE
		BCH_object.decode_bch();
		// print out original and decoded Data
		std::cout<<std::endl<<"Results: "<<std::endl;
		std::cout<<"Original Data  = ";
		for (int i = BCH_object.k-1; i >= 0; i--) {
			std::cout<<BCH_object.Data[i];
		}
		std::cout<<std::endl<<"Recovered Data = ";
		for (int i = BCH_object.n - BCH_object.k; i < BCH_object.n; i++) {
			std::cout<<BCH_object.recD[i];
		}
		// decoding errors: we compare only the Data portion 
		std::cout<<std::endl<<"                 ";
		for (int i = BCH_object.n - BCH_object.k; i < BCH_object.n; i++) {
			if (BCH_object.Data[BCH_object.n - i-1] != BCH_object.recD[i]) {
				BCH_object.decerror++;
				std::cout<<'^';
			} else {
				std::cout<<' ';
			}
		}
		if (BCH_object.decerror) {
			std::cout<<std::endl<<BCH_object.decerror<<" Message decoding errors (at ^)\n\n\n";
			BCH_object.decerror = 0;
		} else {
			std::cout<<std::endl<<"Succesful decoding\n\n\n";
		}
		std::cout<<"Run program again? (y/n)"<<std::endl;
		std::cin>>run_program;
		while(!std::cin || run_program != 'y' && run_program != 'n') {
			std::cout<<"Run program again? (y/n)\r"<<std::endl;
			BCH_object.cin_clean();
			std::cin>>run_program;
		}
	}
    return 0;
}