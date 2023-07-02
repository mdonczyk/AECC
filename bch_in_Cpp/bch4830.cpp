#include "bch4830.hpp"

void bch::read_p() {
	p = 0b1011011;
	primitive_polynomial = p.to_ulong();
	std::cout << "Primitive polynomial:" << std::endl << "p(x) = ";
	verbose_polynomial(p);
}

inline int bch::MSB(const std::bitset <n> &Polynomial) {
	return (n + (GFB-n) - std::countl_zero(Polynomial.to_ullong()));
}

void bch::generate_gf() {
	bch::index_of[0] = -1;
	m = bch::MSB(p);
	for (int i = 0; i < GFB; i++) {
		if (i < m) {
			bch::alpha_to[i] = 1 << i;
			bch::index_of[bch::alpha_to[i]] = i;
		} else {
			if (bch::alpha_to[i - 1] >= 32) {
				bch::alpha_to[i] = (bch::alpha_to[i - 1] << 1) ^ primitive_polynomial;
			} else {
				bch::alpha_to[i] = bch::alpha_to[i - 1] << 1;
			}
			bch::index_of[bch::alpha_to[i]] = i;
		}
	}
}

void bch::gen_poly() {
	std::vector <std::vector <int>> cycle_cosets;
	std::set <int> unique_elements;
	std::pair<std::set<int>::iterator, bool> status;
	int coset_element = 0;
	cycle_cosets.push_back(std::vector<int>{0});
	while (unique_elements.size() < GFB) {
		status.second = false;
		while (!status.second) {
				coset_element++;
				status = unique_elements.emplace(coset_element);
			}
		std::vector <int> coset;
		coset.push_back(coset_element);
		for (int i = 1; i < m; i++) {
			coset.push_back((coset[i-1] << 1) % GFB);
			status = unique_elements.emplace(coset[i]);
			if (!status.second) {
				break;
			}
			if (coset[0] == (coset[i] << 1) % GFB) {
				cycle_cosets.push_back(coset);
				break;
			}
		}
	}
	for (const auto& coset : cycle_cosets) {
		bool root_found = false;
		for (const auto& element : coset) {
			for (int root = 1; root <= 2*t; root++) {
				if (element == root) {
					root_found = true;
					break;
				}
			}
			if(root_found) {
				zeros_cosets.push_back(coset);
				break;
			}
		}
		if (zeros_cosets.size() == t) {
			break;
		}
	}
	std::vector <int> min_polynomials;
	for (const auto& zero_coset : zeros_cosets) {
		int product = bch::alpha_to[zero_coset[0]] ^ 2; // (ax + x)
		for (uint i = 1; i < zero_coset.size(); i++) {
			product = multiply_int_polynomials(product, bch::alpha_to[zero_coset[i]] ^ 2);
		}
		product %= GFB+1;
		product ^= primitive_polynomial;
		min_polynomials.push_back(product);
	}
	int generator_polynomial = min_polynomials[0];
	for (uint i=1; i<t; i++) {
		generator_polynomial = multiply_int_polynomials(generator_polynomial, min_polynomials[i]);
	}
	generator_polynomial_bitset = generator_polynomial;
	reverse_bitset(generator_polynomial_bitset, k-1);
	std::cout << "This is a (" << n << "," << k << "," << t*2+1 << ") binary bch code" << std::endl;
	std::cout << "g(x) is " << generator_polynomial_bitset.to_string().substr(n - (n-k+1)) << std::endl;
}

void BCH_code_short_t3::encode_bch() {
	std::bitset<n> codeword;
	std::bitset<n> shifted_data = this->Data<<(n-k);
	std::bitset<n> rb = divide_bitset_polynomials(shifted_data, bch::generator_polynomial_bitset).first;
	codeword = shifted_data ^ rb;
	if (verbose_flag) {
		std::cout << std::endl;
	}
	this->Codeword = codeword;
}

std::vector <int> BCH_code_short_t3::calculate_syndromes(const std::bitset<n> &Received_Codeword, bool &errors_in_codeword) {
	std::vector <int> syndromes(2*t+1);
	for (int i = 1; i <= 2*t; i++) {
		syndromes[i] = 0;
		for (int j = 0; j < n; j++) {
			if (Received_Codeword[n - j - 1] != 0) {
				syndromes[i] ^= bch::alpha_to[(i*j) % GFB];
			}
		}
		syndromes[i] = bch::index_of[syndromes[i]];
		// when there are no errors all syndromes should be -1
		if (syndromes[i] != -1) {
			errors_in_codeword = true;
		}
	}
	if (verbose_flag) {
		std::cout << "syndromes= ( ";
		for (int i = 1; i <= 2*t; i++) {
			std::cout << syndromes[i] << " ";
		}
		std::cout << ")" << std::endl;
	}
	return syndromes;
}

Status BCH_code_short_t3::decode_bch() {
	/*
	 Simon Rockliff's implementation of Berlekamp's algorithm.
	 Assume we have received bits in Received_Codeword[i], i=0..(n-1).
	
	 Compute the 2*t syndromes by substituting alpha^i into rec(X) and
	 evaluating, storing the syndromes in syndromes[i], i=1..2t (leave syndromes[0] zero) .
	 Then we use the Berlekamp algorithm to find the error location polynomial
	 elp[i].
	
	 If the degree of the elp is >t, then we cannot correct all the errors, and
	 we have detected an uncorrectable error pattern. We output the information
	 bits uncorrected.
	
	 If the degree of elp is <=t, we substitute alpha^i , i=1..n into the elp
	 to get the roots, hence the inverse roots, the error location numbers.
	 This step is usually called "Chien's search".
	
	 If the number of errors located is not equal the degree of the elp, then
	 the decoder assumes that there are more than t errors and cannot correct
	 them, only detect them. We output the information bits uncorrected.
	*/

	bool errors_in_codeword = false;
	std::vector <int> syndromes = calculate_syndromes(this->Received_Codeword, errors_in_codeword);
	std::bitset <n> Decoded_Codeword = Received_Codeword;
	if (errors_in_codeword) {
		/*
		 Compute the error location polynomial via the Berlekamp
		 iterative algorithm. Following the terminology of Lin and
		 Costello's book :   d[u] is the 'mu'th discrepancy, where
		 u='mu'+1 and 'mu' is the step number
		 ranging from -1 to 2*t (see L&C),  l[u] is the degree of
		 the elp at that step, and u_l[u] is the difference between
		 the step number and the degree of the elp. 

		*/
		std::vector <int> error_locations;
		int elp[100][100] = {{0, 0}}; // error location polynomial
		int d[100] = {0}, l[100] = {0}, u_lu[100] = {0};
		for (int i = 1; i < t*2; i++) {
			elp[0][i] = -1; // index form
			elp[1][i] = 0; // polynomial form
		}
		u_lu[0] = -1;
		u_lu[1] = 0;
		d[1] = syndromes[1]; // index form
		elp[1][0] = 1; // polynomial form
		int u=0;
		do {
			u++;
			if (d[u] == -1) {
				l[u + 1] = l[u];
				for (int i = 0; i <= l[u]; i++) {
					elp[u + 1][i] = elp[u][i];
					elp[u][i] = bch::index_of[elp[u][i]];
				}
			} else {
				// search for words with greatest u_lu[q] for which d[q]!=0 
				int q = u - 1;
				while ((d[q] == -1) && (q > 0)) {
					q--;
				}
				// have found first non-zero d[q]
				if (q > 0) {
				  int j = q;
				  do {
				    j--;
				    if ((d[j] != -1) && (u_lu[q] < u_lu[j]))
				      q = j;
				  } while (j > 0);
				}
				// have now found q such that d[u]!=0 and
				// u_lu[q] is maximum 
				// store degree of new elp polynomial
				if (l[u] > l[q] + u - q) {
					l[u + 1] = l[u];
				} else {
					l[u + 1] = l[q] + u - q;
				}
				// form new elp(x)
				for (int i = 0; i < 2*t; i++) {
					elp[u + 1][i] = 0;
				}
				for (int i = 0; i <= l[q]; i++) {
					if (elp[q][i] != -1) {
						elp[u + 1][i + u - q] = bch::alpha_to[(d[u] + GFB - d[q] + elp[q][i]) % GFB];
					}
				}
				for (int i = 0; i <= l[u]; i++) {
					elp[u + 1][i] ^= elp[u][i];
					elp[u][i] = bch::index_of[elp[u][i]];
				}
			}
			u_lu[u + 1] = u - l[u + 1];
			// form (u+1)th discrepancy
			if (u < 2*t) {
			// no discrepancy computed on last iteration
				if (syndromes[u + 1] != -1) {
					d[u + 1] = bch::alpha_to[syndromes[u + 1]];
				} else {
				d[u + 1] = 0;
				}
				for (int i = 1; i <= l[u + 1]; i++) {
					if ((syndromes[u + 1 - i] != -1) && (elp[u + 1][i] != 0)) {
						d[u + 1] ^= bch::alpha_to[(syndromes[u + 1 - i] + bch::index_of[elp[u + 1][i]]) % GFB];
					}
				}
				// put d[u+1] into index form
				d[u + 1] = bch::index_of[d[u + 1]];
			}
		} while ((u < 2*t) && (l[u + 1] <= t));
		u++;
		if (l[u] <= t) {// Can correct errors
			// put elp into index form
			if (verbose_flag) {
				std::cout << "Sigma(x) = ";
			}
			for (int i = 0; i <= l[u]; i++) {
				elp[u][i] = bch::index_of[elp[u][i]];
				if (verbose_flag) {
					std::cout << elp[u][i] << " ";
				}
			}
			// Find roots of the error location polynomial:
			// Chien search
			if (verbose_flag) {
				std::cout << std::endl << "Roots: ";
			}
			for (int i = 1; i <= GFB; i++) {
				int q = 1;
				for (int j = 1; j <= l[u]; j++) {
					elp[u][j] = (elp[u][j] + j) % GFB;
					q ^= bch::alpha_to[elp[u][j]];
				}
				// Store error location number indices
				if (!q) {
					error_locations.push_back(n - i + (GFB-n));
					if (verbose_flag) {
						std::cout << error_locations.back() << " ";
					}
				}
			}
			if (verbose_flag) {
				std::cout << std::endl;
			}
			if (error_locations.size() == (unsigned)l[u]) {
			// no. roots = degree of elp hence <= t errors
				for (auto const &error_location : error_locations) {
					auto err_loc = (n-1-error_location + n) % n;
					if (err_loc < 0 || err_loc >= n) {
						if (verbose_flag) {
							std::cout<<std::endl<<"Incomplete decoding: errors detected (err_loc out of bound: err_loc = " << err_loc << ")"<<std::endl;
						}
						this->Decoded_Data = (std::bitset<k>)Decoded_Codeword.to_string().substr(0, k);
						return FAIL;
					} else {
						Decoded_Codeword.flip(err_loc);
					}
				}
			} else {// elp has degree > t hence cannot solve
				if (verbose_flag) {
					std::cout<<std::endl<<"Incomplete decoding: errors detected (elp has degree > t) "<<std::endl;
				}
				this->Decoded_Data = (std::bitset<k>)Decoded_Codeword.to_string().substr(0, k);
				return FAIL;
			}
		} else {
			if (verbose_flag) {
				std::cout<<std::endl<<"Incomplete decoding: errors detected (l[u] > t: l[u] = " << l[u] << ")"<<std::endl;
			}
			this->Decoded_Data = (std::bitset<k>)Decoded_Codeword.to_string().substr(0, k);
			return FAIL;
		}
	} else {
		if (verbose_flag) {
			std::cout << "No errors found" << std::endl;
		}
	}
	this->Decoded_Data = (std::bitset<k>)Decoded_Codeword.to_string().substr(0, k);
	return SUCCESS;
}

void bch::verbose_polynomial(const std::bitset <n> &Polynomial) {
	int power = bch::MSB(Polynomial);
	for (int i=power; i>=0; i--) {
		if (Polynomial[i]) {
			if (i != power) {
				std::cout << " + " ;
			}
			if (i!=0) {
				if (i == 1) {
					std::cout << "x";
				} else {
				std::cout << "x^" << i;
				}
			} else {
				std::cout << "1";
			}
		}
	}
	std::cout << std::endl;
}

template<size_t N>
void bch::reverse_bitset(std::bitset <N> &Polynomial, int Shift) {
    for(size_t i = 0; i < N/2; ++i) {
    	bool temp_bit = Polynomial[i];
    	Polynomial[i] = Polynomial[N-i-1];
    	Polynomial[N-i-1] = temp_bit;
    }
	Polynomial >>= Shift;
}

std::pair<std::bitset <n>, std::bitset <n>> BCH_code_short_t3::divide_bitset_polynomials(
			const std::bitset <n> &Dividend, const std::bitset <n> &Divisor) {
	std::bitset <n> quotient, remainder = Dividend;
 	while (bch::MSB(remainder) >= bch::MSB(Divisor)) {
		int shift = bch::MSB(remainder) - bch::MSB(Divisor);
		remainder ^= Divisor << shift;
		quotient.flip(shift); 
	}
	return {remainder, quotient};
}

int bch::multiply_int_polynomials(int Mulitplicand, int Multiplicator) {
	int product = 0;
	while (Mulitplicand > 0) {
		if (Mulitplicand & 1) {
			product ^= Multiplicator;
		}
		Multiplicator <<= 1;
		Mulitplicand >>= 1;
	}
	return product;
}

void BCH_code_short_t3::print_original_codeword_and_received_codeword() {
	if (verbose_flag) {
		std::cout << "c(x) = "<< this->Codeword << std::endl;
		std::cout << "r(x) = "<< this->Received_Codeword << std::endl;
	}
	std::bitset <n> bit_difference = Codeword ^ Received_Codeword;
	if (bit_difference.count() != 0) {
		if (verbose_flag) {
			std::cout<<"       "<<(bit_difference.to_string(' ','^'))<<std::endl;
			if (bit_difference.count() == 1) {
				std::cout << "Position of " << bit_difference.count() <<
				 " error in the received codeword at ^." << std::endl;
			} else {
				std::cout << "Positions of " << bit_difference.count() <<
				 " errors in the received codeword at ^." << std::endl;
			}
		}
		if (bit_difference.count() > t) {
			g_big_errors_count++;
		}
	} else {
		if (verbose_flag) {
			std::cout << "No errors in the received codeword." << std::endl;
		}
	}
}

void BCH_code_short_t3::print_original_message_and_decoded_message() {
	if (verbose_flag) {
		std::cout <<  DASH_LINE << "Results" << DASH_LINE << std::endl;
	}
	std::bitset<k> Original_Data(this->Data.to_string().substr(n-k, n));
	if (verbose_flag) {
		std::cout << "Original Data  = " << Original_Data << std::endl;
		std::cout << "Recovered Data = " << this->Decoded_Data<< std::endl;
	}
	std::bitset <k> bit_difference = Original_Data ^ Decoded_Data;
	if (bit_difference.count()) {
		if (verbose_flag) {
			std::cout << "                 " << (bit_difference.to_string(' ','^')) << std::endl;
			if (bit_difference.count() == 1) {
				std::cout << "Position of " << bit_difference.count() << " message decoding error at ^." << std::endl;
			} else {
				std::cout << "Positions of " << bit_difference.count() << " message decoding errors at ^." << std::endl;
			}
		}
		g_uncaught_errors_count++;
	} else {
		if (verbose_flag) {
			std::cout <<"Successful decoding." << std::endl;
		}
	}
}

void BCH_code_short_t3::introduce_errors() {
	std::bitset <n> received_codeword = this->Codeword;
	std::random_device rd;
	std::mt19937 gen (rd());
	std::uniform_int_distribution<> d (0, bch::error_probability);
	for (int i = n - 1; i >= 0; i--) {
		int randnum = d(gen);
		if (randnum == 0) {
			g_introduced_errors_count ++;
			received_codeword.flip(i);
		}
	}
	this->Received_Codeword = received_codeword;
}

void print_help_message(const char *file_name) {
	std::cout << "Usage:\n"
		<< file_name << " [-h] -i image -p err_prob [-v]\n\n"
		<< "Options:\n"
		<< "  -i image			Choose one image from images folder, example: images/image2.bmp.\n"
		<< "  -p err_prob		Give probability between (10 and 10000000) that a 1 in err_prob error will occur in \n"
		<< "			   the codeword during a simulated transmission through a noisy medium, example: 1000.\n"
		<< "Optional arguments:\n"
		<< "  -h	Show this help message.\n"
		<< "  -v	Enable verbose_flag encoding and decoding logs which will print out the whole process to the \n"
		<< "	   terminal, is disabled by default. WARNING! This option causes the threads to run sequentially instead \n"
		<< " 	   of in parallel which combined with printing operations to console causes a severe performance degradation.\n";
}

void begin_main_process(const int thread_id, const int MBTG_beginning, const int MBTG_end) {
	if (verbose_flag) { bch::Mutex.lock(); }
	for (int i = MBTG_beginning; i < MBTG_end; i++) {
		if (verbose_flag) {
			std::cout << LINE << " Worker " << thread_id << " START " << LINE;
		}
		
		bch::BCH_objects[i] = make_shared<BCH_code_short_t3>(std::bitset<n>(bch::vector_of_message_polynomials[i].to_string()));

		bch::BCH_objects[i]->encode_bch();

		bch::BCH_objects[i]->introduce_errors();

		bch::BCH_objects[i]->print_original_codeword_and_received_codeword();

		bch::BCH_objects[i]->Received_Data = std::bitset <k>(bch::BCH_objects[i]->Received_Codeword.to_string().substr(0, k));

		bool status = bch::BCH_objects[i]->decode_bch();

		// check the decoding status
		if (status == SUCCESS) {
			bch::BCH_objects[i]->print_original_message_and_decoded_message();
			g_success_count++;
		} else {
			g_failure_count++;
		}

		if (verbose_flag) {
			std::cout << LINE << " Worker " << thread_id << " STOP *" << LINE << std::endl << std::endl;
		}
	}
	if (verbose_flag) { bch::Mutex.unlock(); }
}

void populate_vector_of_message_polynomials (const unsigned char* str, const int MBTG_beginning, const int MBTG_end, 
											const int MPTG_end, const int Bit_pos) {
	std::bitset<k> temp_bitset;
	int it = MPTG_end;
	int bit_pos = Bit_pos;
	for (int i = MBTG_beginning; i < MBTG_end; i++) {
		for (int j = 0; j < 8; ++j) {
			temp_bitset[bit_pos] = str[i] & (1 << j);
			bit_pos++;
			if (bit_pos == k) {
				bch::vector_of_message_polynomials[it] ^= temp_bitset;
				it++;
				temp_bitset.reset();
				bit_pos = 0;
			}
		}
	}
	if (bit_pos != 0) {
		bch::vector_of_message_polynomials[it] ^= temp_bitset;
	}
}

void populate_unsigned_char_vectors (const int MBTG_beginning, const int MBTG_end, 
										const int MPTG_end, const int Bit_pos) {
	int it = MPTG_end;
	int bit_pos = Bit_pos;
	for (int i = MBTG_beginning; i < MBTG_end; i++) {
		for (int j = 0; j < 8; ++j) {
			bch::recovered_charstream[i] ^= bch::BCH_objects[it]->Decoded_Data[bit_pos] << j;
			bch::modified_charstream[i] ^= bch::BCH_objects[it]->Received_Data[bit_pos] << j;
			bit_pos++;
			if (bit_pos == k) {
				it++;
				bit_pos = 0;
			}
		}
	}
}

int parse_arguments (const int argc, char* argv[]) {

	int fd = FAIL;
	int c;

	if (argc != 2 && argc != 5 && argc != 6) {
		std::cout << "Invalid number of arguments"<< std::endl;
		return FAIL;
	}

	while ((c = getopt (argc, argv, "i:p:vh")) != FAIL) {
		switch (c) {
			case 'h':
				print_help_message(argv[0]);
				exit(0);
			case 'i':
				bch::filename = optarg;
				fd = open(optarg, O_RDONLY);
				if (fd == FAIL) {
					std::cout << "Incorrect Image file name" << std::endl;
					return FAIL;
				}
				break;
			case 'p':
				bch::error_probability = atoi(optarg);
				if (bch::error_probability < 10 || bch::error_probability > 10000000) {
					std::cout << "Error probability should be between (10 and 10000000)" << std::endl;
					return FAIL;
				}
				break;
			case 'v':
				verbose_flag = true;
				break;
			case '?':
					std::cout << optopt << std::endl;
				return FAIL;
			default:
				abort();
		}
	}
	
	return fd;
}

int main(int argc, char* argv[]) {
	int fd {parse_arguments(argc, argv)};
	
	if (fd == FAIL) {
		std::cout << "Argument parsing failed, type: \"" << argv[0] << " -h\" for more information" << std::endl;
		return 1;
	}
	
	struct stat sb;

	if (fstat(fd, &sb) == FAIL) {
		std::cout << "Failed to get file size\n";
		return 1;
	}

	const int FILE_BYTE_SIZE = sb.st_size;
	
	unsigned char* buffer = (unsigned char*) mmap(
		NULL, 
		FILE_BYTE_SIZE, 
		PROT_READ, MAP_PRIVATE, 
		fd, 
		0);

	std::string image_with_errors_path = bch::filename.substr(0, bch::filename.length()-4).append("_with_errors_BCH_code_short_t3.bmp");
	std::string image_fixed_path = bch::filename.substr(0, bch::filename.length()-4).append("_fixed_BCH_code_short_t3.bmp");

	remove(image_with_errors_path.c_str());
	remove(image_fixed_path.c_str());

	std::ofstream image_with_errors (image_with_errors_path.c_str(), std::ios::out | std::ios::app | std::ios::binary);
	std::ofstream image_fixed (image_fixed_path.c_str(), std::ios::out | std::ios::app | std::ios::binary);

	// do a check to see if there are any std::left over bites and if there are we need to add another message polynomial
	const ssize_t NUMBER_OF_MESSAGE_POLYNOMIALS = (FILE_BYTE_SIZE*8)%k? 
													(FILE_BYTE_SIZE*8)/k + 1 :
													(FILE_BYTE_SIZE*8)/k;

	bch::vector_of_message_polynomials.resize(NUMBER_OF_MESSAGE_POLYNOMIALS);
	bch::BCH_objects.resize(NUMBER_OF_MESSAGE_POLYNOMIALS);
	bch::recovered_charstream.resize(FILE_BYTE_SIZE);
	bch::modified_charstream.resize(FILE_BYTE_SIZE);

	std::vector<std::thread> threads;
	const int NUM_THREADS = std::thread::hardware_concurrency();
	std::cout << "number of detected threads: " << NUM_THREADS << std::endl;

	// MESSAGE_POLYNOMIALS_THREAD_GROUP --> count of group of message polynomials that will be used by a single thread
	const int MESSAGE_POLYNOMIALS_THREAD_GROUP = NUMBER_OF_MESSAGE_POLYNOMIALS%NUM_THREADS? 
													NUMBER_OF_MESSAGE_POLYNOMIALS/NUM_THREADS :
													NUMBER_OF_MESSAGE_POLYNOMIALS/NUM_THREADS - 1;

	// MESSAGE_BYTES_THREAD_GROUP --> count of group of message bytes that will be used by a single thread
	const int MESSAGE_BYTES_THREAD_GROUP = FILE_BYTE_SIZE / NUM_THREADS;

	std::cout << "Parsing image file..." << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	int bit_pos;
	int old_bit_pos{0};
	int overlaping_message_polynomial{0};

	// MESSAGE_BYTES_THREAD_GROUP beginning
	int MBTG_beginning;
	// MESSAGE_BYTES_THREAD_GROUP end
	int MBTG_end;
	// MESSAGE_POLYNOMIALS_THREAD_GROUP end
	int MPTG_end;

	for (int thread_id = 0; thread_id < NUM_THREADS; thread_id++) {
		MBTG_beginning = MESSAGE_BYTES_THREAD_GROUP*thread_id;
		MBTG_end = MESSAGE_BYTES_THREAD_GROUP*(thread_id + 1);

		if (thread_id == NUM_THREADS-1) {
			MBTG_end = FILE_BYTE_SIZE;
		}

		bit_pos = (MBTG_beginning*8) % k;

		if (bit_pos < old_bit_pos) { overlaping_message_polynomial++; }

		MPTG_end = MESSAGE_POLYNOMIALS_THREAD_GROUP*thread_id + overlaping_message_polynomial;

		old_bit_pos = bit_pos;

		threads.push_back(std::thread(populate_vector_of_message_polynomials, buffer, MBTG_beginning, 
									MBTG_end, MPTG_end, bit_pos));
	}

	for (auto& thread : threads) {
		thread.join();
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto dur = stop - start;
	auto duration = std::chrono::duration_cast<std::chrono::duration<float>>(dur);

	std::cout << "Parsing image file done and it took: ";
	std::cout << std::to_string(duration.count()).substr(0, std::to_string(duration.count()).length()-3) << " seconds" << std::endl;

	start = std::chrono::high_resolution_clock::now();

	bch::initialize_BCH();

	threads.clear();

	if (!verbose_flag){
		std::cout<<"Please be patient, starting coding and decoding process..." << std::endl;
	}

	for (int thread_id = 0; thread_id < NUM_THREADS; thread_id++) {
		MBTG_beginning = MESSAGE_POLYNOMIALS_THREAD_GROUP*thread_id;
		MBTG_end = MESSAGE_POLYNOMIALS_THREAD_GROUP*(thread_id + 1);

		if (thread_id == NUM_THREADS-1) {
			MBTG_end = NUMBER_OF_MESSAGE_POLYNOMIALS;
		}

		threads.push_back(std::thread(begin_main_process, thread_id, MBTG_beginning, MBTG_end));
	}

	for (auto& thread : threads) {
		thread.join();
	}

	stop = std::chrono::high_resolution_clock::now();
	dur = stop - start;
	auto main_duration = std::chrono::duration_cast<std::chrono::duration<float>>(dur);

	std::cout << "Coding and decoding process done and it took: ";
	std::cout << std::to_string(main_duration.count()).substr(0, std::to_string(main_duration.count()).length()-3) << " seconds" << std::endl;

	std::cout << "Converting modified and recovered data from bitsets to bytes..." << std::endl;

	start = std::chrono::high_resolution_clock::now();

	threads.clear();

	old_bit_pos = 0;
	overlaping_message_polynomial = 0;

	for (int thread_id = 0; thread_id < NUM_THREADS; thread_id++) {
		MBTG_beginning = MESSAGE_BYTES_THREAD_GROUP*thread_id;
		MBTG_end = MESSAGE_BYTES_THREAD_GROUP*(thread_id + 1);

		if (thread_id == NUM_THREADS-1) {
			MBTG_end = FILE_BYTE_SIZE;
		}

		bit_pos = (MBTG_beginning*8) % k;

		if (bit_pos < old_bit_pos) { overlaping_message_polynomial++; }

		MPTG_end = MESSAGE_POLYNOMIALS_THREAD_GROUP*thread_id + overlaping_message_polynomial;

		old_bit_pos = bit_pos;

		threads.push_back(std::thread(populate_unsigned_char_vectors, MBTG_beginning, MBTG_end, 
									MPTG_end, bit_pos));
	}

	for (auto& thread : threads) {
		thread.join();
	}

	stop = std::chrono::high_resolution_clock::now();
	dur = stop - start;
	duration = std::chrono::duration_cast<std::chrono::duration<float>>(dur);

	std::cout << "Converting modified and recovered data from bitsets to bytes done and it took: ";
	std::cout << std::to_string(duration.count()).substr(0, std::to_string(duration.count()).length()-3) << " seconds" << std::endl;

	for (int i = 0; i < HEADER_BYTES; i++) {
		image_with_errors << buffer[i];
		image_fixed << buffer[i];
	}

	for (int i = HEADER_BYTES; i < FILE_BYTE_SIZE; i++) {
		image_with_errors << bch::modified_charstream[i];
		image_fixed << bch::recovered_charstream[i];
	}

	int difference_count = 0;
	for (int i = 0; i < NUMBER_OF_MESSAGE_POLYNOMIALS; i++) {
		difference_count += (bch::vector_of_message_polynomials[i] ^ bch::BCH_objects[i]->Decoded_Data).count();
	}

	// int total_fresh_diff = accumulate (bch::BCH_objects.begin(), bch::BCH_objects.end(), 0, 
	//  	[](int sum, shared_ptr<BCH_code_short_t3> const& obj) 
	// 	{ return sum + (obj->Codeword ^ obj->Received_Codeword).count(); });
	// if (total_fresh_diff != g_introduced_errors_count) { std::cout<< "Some of introduced errors cancelled each other out" << std::endl; }

	auto table_row_separator = []()->void { std::cout << std::setw(42) << std::setfill('-') << "+" << std::setw(20) << "+" << std::setw(1) << "+" << std::setfill (' ') << std::endl; };
	auto table_items = [](const std::string &lhs, const std::string &rhs)->void { std::cout << "| " << std::setw(40) << lhs << std::setw(20) << "| " + rhs << "|" << std::endl; };

	std::cout << std::endl << std::left;
	table_row_separator();
	table_items("Code used ", "(" + std::to_string(n) + "," + std::to_string(k) + "," + std::to_string(t*2+1 ) + ")");
	table_row_separator();
	table_items("Original image used ", bch::filename);
	table_row_separator();
	table_items("Given error probability ", std::to_string(1/(float)bch::error_probability));
	table_row_separator();
	table_items("Real error probability ", std::to_string(1/(float)((NUMBER_OF_MESSAGE_POLYNOMIALS*n)/g_introduced_errors_count)));
	table_row_separator();
	table_items("Number of all data bits ", std::to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * k));
	table_row_separator();
	table_items("Number of all data + redundant bits ", std::to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * n));
	table_row_separator();
	table_items("Number of all introduced errors ", std::to_string(g_introduced_errors_count));
	table_row_separator();
	table_items("Number of all message polynomials ", std::to_string(NUMBER_OF_MESSAGE_POLYNOMIALS));
	table_row_separator();
	table_items("Number of successful decodings", std::to_string(g_success_count));
	table_row_separator();
	table_items("Number of decoding errors ", std::to_string(g_failure_count));
	table_row_separator();
	table_items("Number of uncaught decoding errors ", std::to_string(g_uncaught_errors_count));
	table_row_separator();
	table_items("Number of over t errors in codewords ", std::to_string(g_big_errors_count));
	table_row_separator();
	table_items("Final data bits difference ", std::to_string(difference_count));
	table_row_separator();
	table_items("Encoding and decoding time ", std::to_string(main_duration.count()).substr(0, std::to_string(main_duration.count()).length()-3) + " seconds ");
	table_row_separator();

	std::cout << std::endl << "Results have been written to BCH_logs.txt file" << std::endl;
	std::cout << std::endl << "To view made images on linux use \"feh -F -Z --force-aliasing -d " << bch::filename << " " << image_with_errors_path.c_str() << " " << image_fixed_path.c_str() << "\"" << std::endl << std::endl;

	std::ofstream BCH_logs;
	BCH_logs.open ("BCH_logs.txt", std::ios::out | std::ios::app);
	auto BCH_logs_items = [&](const std::string &lhs, const std::string &rhs)->void { BCH_logs << std::setw(37) << lhs << std::setw(20) << "| " + rhs << std::endl; };

	BCH_logs << std::endl << std::left;
	BCH_logs_items("Code used ", "(" + std::to_string(n) + "," + std::to_string(k) + "," + std::to_string(t*2+1 ) + ")");
	BCH_logs_items("Original image used ", bch::filename);
	BCH_logs_items("Given error probability ", std::to_string(1/(float)bch::error_probability));
	BCH_logs_items("Real error probability ", std::to_string(1/(float)((NUMBER_OF_MESSAGE_POLYNOMIALS*n)/g_introduced_errors_count)));
	BCH_logs_items("Number of all data bits ", std::to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * k));
	BCH_logs_items("Number of all data + redundant bits ", std::to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * n));
	BCH_logs_items("Number of all generated errors ", std::to_string(g_introduced_errors_count));
	BCH_logs_items("Number of all message polynomials ", std::to_string(NUMBER_OF_MESSAGE_POLYNOMIALS));
	BCH_logs_items("Number of successful decodings ", std::to_string(g_success_count));
	BCH_logs_items("Number of decoding errors ", std::to_string(g_failure_count));
	BCH_logs_items("Number of uncaught decoding errors", std::to_string(g_uncaught_errors_count));
	BCH_logs_items("Number of over t errors in codewords", std::to_string(g_big_errors_count));
	BCH_logs_items("Final data bits difference ", std::to_string(difference_count));
	BCH_logs_items("Encoding and decoding time ", std::to_string(main_duration.count()).substr(0, std::to_string(main_duration.count()).length()-3) + " seconds ");

	munmap(buffer, FILE_BYTE_SIZE);
	close(fd);

	image_with_errors.close();
	image_fixed.close();
	BCH_logs.close();

	return 0;
}
