#include "bch4836.hpp"

using namespace std;
using namespace BCH;

void BCH::read_p() {
	p = 0b1011011;
	primitive_polynomial = p.to_ulong();
	cout << "Primitive polynomial:" << endl << "p(x) = ";
	verbose_polynomial(p);
}

int BCH::MSB(const bitset <n> &Polynomial) {
	return (n + (GFB-n) - countl_zero(Polynomial.to_ullong()));
}

void BCH::generate_gf() {
	index_of[0] = -1;
	m = MSB(p);
	for (int i = 0; i < GFB; i++) {
		if (i < m) {
			alpha_to[i] = 1 << i;
			index_of[alpha_to[i]] = i;
		} else {
			if (alpha_to[i - 1] >= 32) {
				alpha_to[i] = (alpha_to[i - 1] << 1) ^ primitive_polynomial;
			} else {
				alpha_to[i] = alpha_to[i - 1] << 1;
			}
			index_of[alpha_to[i]] = i;
		}
	}
}

void BCH::gen_poly() {
	vector <vector <int>> cycle_cosets;
	set <int> unique_elements;
	pair<set<int>::iterator, bool> status;
	int coset_element = 0;
	cycle_cosets.push_back(vector<int>{0});
	while (unique_elements.size() < GFB) {
		status.second = false;
		while (!status.second) {
				coset_element++;
				status = unique_elements.emplace(coset_element);
			}
		vector <int> coset;
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
	vector <int> min_polynomials;
	for (const auto& zero_coset : zeros_cosets) {
		int product = alpha_to[zero_coset[0]] ^ 2; // (ax + x)
		for (uint i = 1; i < zero_coset.size(); i++) {
			product = multiply_int_polynomials(product, alpha_to[zero_coset[i]] ^ 2);
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
	cout << "This is a (" << n << "," << k << "," << t*2+1 << ") binary BCH code" << endl;
	cout << "g(x) is " << generator_polynomial_bitset.to_string().substr(n - (n-k+1)) << endl;
}

void BCH_code_short_t2::encode_bch() {
	bitset<n> codeword;
	bitset<n> shifted_data = this->Data<<(n-k);
	bitset<n> rb = divide_bitset_polynomials(shifted_data, generator_polynomial_bitset).first;
	codeword = shifted_data ^ rb;
	if (verbose_flag) {
		cout << endl;
	}
	this->Codeword = codeword;
}

vector <int> BCH_code_short_t2::calculate_syndromes(const bitset<n> &Received_Codeword, bool &errors_in_codeword) {
	vector <int> syndromes(2*t+1);
	for (int i = 1; i <= 2*t; i++) {
		syndromes[i] = 0;
		for (int j = 0; j < n; j++) {
			if (Received_Codeword[n - j - 1] != 0) {
				syndromes[i] ^= alpha_to[(i*j) % GFB];
			}
		}
		syndromes[i] = index_of[syndromes[i]];
		// when there are no errors all syndromes should be -1
		if (syndromes[i] != -1) {
			errors_in_codeword = true;
		}
	}
	if (verbose_flag) {
		cout << "syndromes= ( ";
		for (int i = 1; i <= 2*t; i++) {
			cout << syndromes[i] << " ";
		}
		cout << ")" << endl;
	}
	return syndromes;
}

Status BCH_code_short_t2::decode_bch() {
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
	vector <int> syndromes = calculate_syndromes(this->Received_Codeword, errors_in_codeword);
	bitset <n> Decoded_Codeword = Received_Codeword;
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
		vector <int> error_locations;
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
					elp[u][i] = index_of[elp[u][i]];
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
						elp[u + 1][i + u - q] = alpha_to[(d[u] + GFB - d[q] + elp[q][i]) % GFB];
					}
				}
				for (int i = 0; i <= l[u]; i++) {
					elp[u + 1][i] ^= elp[u][i];
					elp[u][i] = index_of[elp[u][i]];
				}
			}
			u_lu[u + 1] = u - l[u + 1];
			// form (u+1)th discrepancy
			if (u < 2*t) {
			// no discrepancy computed on last iteration
				if (syndromes[u + 1] != -1) {
					d[u + 1] = alpha_to[syndromes[u + 1]];
				} else {
				d[u + 1] = 0;
				}
				for (int i = 1; i <= l[u + 1]; i++) {
					if ((syndromes[u + 1 - i] != -1) && (elp[u + 1][i] != 0)) {
						d[u + 1] ^= alpha_to[(syndromes[u + 1 - i] + index_of[elp[u + 1][i]]) % GFB];
					}
				}
				// put d[u+1] into index form
				d[u + 1] = index_of[d[u + 1]];
			}
		} while ((u < 2*t) && (l[u + 1] <= t));
		u++;
		if (l[u] <= t) {// Can correct errors
			// put elp into index form
			if (verbose_flag) {
				cout << "Sigma(x) = ";
			}
			for (int i = 0; i <= l[u]; i++) {
				elp[u][i] = index_of[elp[u][i]];
				if (verbose_flag) {
					cout << elp[u][i] << " ";
				}
			}
			// Find roots of the error location polynomial:
			// Chien search
			if (verbose_flag) {
				cout << endl << "Roots: ";
			}
			for (int i = 1; i <= GFB; i++) {
				int q = 1;
				for (int j = 1; j <= l[u]; j++) {
					elp[u][j] = (elp[u][j] + j) % GFB;
					q ^= alpha_to[elp[u][j]];
				}
				// Store error location number indices
				if (!q) {
					error_locations.push_back(n - i + (GFB-n));
					if (verbose_flag) {
						cout << error_locations.back() << " ";
					}
				}
			}
			if (verbose_flag) {
				cout << endl;
			}
			if (error_locations.size() == (unsigned)l[u]) {
			// no. roots = degree of elp hence <= t errors
				for (auto const &error_location : error_locations) {
					auto err_loc = (n-1-error_location + n) % n;
					if (err_loc < 0 || err_loc >= n) {
						if (verbose_flag) {
							cout<<endl<<"Incomplete decoding: errors detected (err_loc out of bound: err_loc = " << err_loc << ")"<<endl;
						}
						this->Decoded_Data = (bitset<k>)Decoded_Codeword.to_string().substr(0, k);
						return FAIL;
					} else {
						Decoded_Codeword.flip(err_loc);
					}
				}
			} else {// elp has degree > t hence cannot solve
				if (verbose_flag) {
					cout<<endl<<"Incomplete decoding: errors detected (elp has degree > t) "<<endl;
				}
				this->Decoded_Data = (bitset<k>)Decoded_Codeword.to_string().substr(0, k);
				return FAIL;
			}
		} else {
			if (verbose_flag) {
				cout<<endl<<"Incomplete decoding: errors detected (l[u] > t: l[u] = " << l[u] << ")"<<endl;
			}
			this->Decoded_Data = (bitset<k>)Decoded_Codeword.to_string().substr(0, k);
			return FAIL;
		}
	} else {
		if (verbose_flag) {
			cout << "No errors found" << endl;
		}
	}
	this->Decoded_Data = (bitset<k>)Decoded_Codeword.to_string().substr(0, k);
	return SUCCESS;
}

void BCH::verbose_polynomial(const bitset <n> &Polynomial) {
	int power = MSB(Polynomial);
	for (int i=power; i>=0; i--) {
		if (Polynomial[i]) {
			if (i != power) {
				cout << " + " ;
			}
			if (i!=0) {
				if (i == 1) {
					cout << "x";
				} else {
				cout << "x^" << i;
				}
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

pair<bitset <n>, bitset <n>> BCH_code_short_t2::divide_bitset_polynomials(
			const bitset <n> &Dividend, const bitset <n> &Divisor) {
	bitset <n> quotient, remainder = Dividend;
 	while (MSB(remainder) >= MSB(Divisor)) {
		int shift = MSB(remainder) - MSB(Divisor);
		remainder ^= Divisor << shift;
		quotient.flip(shift); 
	}
	return {remainder, quotient};
}

int BCH::multiply_int_polynomials(int Mulitplicand, int Multiplicator) {
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

void BCH_code_short_t2::print_original_codeword_and_received_codeword() {
	if (verbose_flag) {
		cout << "c(x) = "<< this->Codeword << endl;
		cout << "r(x) = "<< this->Received_Codeword << endl;
	}
	bitset <n> bit_difference = Codeword ^ Received_Codeword;
	if (bit_difference.count() != 0) {
		if (verbose_flag) {
			cout<<"       "<<(bit_difference.to_string(' ','^'))<<endl;
			if (bit_difference.count() == 1) {
				cout << "Position of " << bit_difference.count() <<
				 " error in the received codeword at ^." << endl;
			} else {
				cout << "Positions of " << bit_difference.count() <<
				 " errors in the received codeword at ^." << endl;
			}
		}
		if (bit_difference.count() > t) {
			g_big_errors_count++;
		}
	} else {
		if (verbose_flag) {
			cout << "No errors in the received codeword." << endl;
		}
	}
}

void BCH_code_short_t2::print_original_message_and_decoded_message() {
	if (verbose_flag) {
		cout <<  DASH_LINE << "Results" << DASH_LINE << endl;
	}
	bitset<k> Original_Data(this->Data.to_string().substr(n-k, n));
	if (verbose_flag) {
		cout << "Original Data  = " << Original_Data << endl;
		cout << "Recovered Data = " << this->Decoded_Data<< endl;
	}
	bitset <k> bit_difference = Original_Data ^ Decoded_Data;
	if (bit_difference.count()) {
		if (verbose_flag) {
			cout << "                 " << (bit_difference.to_string(' ','^')) << endl;
			if (bit_difference.count() == 1) {
				cout << "Position of " << bit_difference.count() << " message decoding error at ^." << endl;
			} else {
				cout << "Positions of " << bit_difference.count() << " message decoding errors at ^." << endl;
			}
		}
		g_uncaught_errors_count++;
	} else {
		if (verbose_flag) {
			cout <<"Successful decoding." << endl;
		}
	}
}

void BCH_code_short_t2::introduce_errors() {
	bitset <n> received_codeword = this->Codeword;
	random_device rd;
	mt19937 gen (rd());
	uniform_int_distribution<> d (0, error_probability);
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
	cout << "Usage:\n"
		<< file_name << " [-h] -i <image> -p err_prob [-v]\n\n"
		<< "Options:\n"
		<< "  -i <file_name>	Choose one image from images folder, example: images/image2.bmp.\n"
		<< "  -p err_prob		Give probability between (10 and 10000000) that a 1 in err_prob error will occur in \n"
		<< "			   the codeword during a simulated transmission through a noisy medium, example: 1000.\n"
		<< "Optional arguments:\n"
		<< "  -h	Show this help message.\n"
		<< "  -v	Enable verbose_flag encoding and decoding logs which will print out the whole process to the \n"
		<< "	   terminal, is disabled by default. WARNING! This option causes the threads to run sequenitally instead \n"
		<< " 	   of in parallel which combined with printing operations to console causes a severe performance degradation.\n";
}

void begin_main_process(const int thread_id, const int thread_zone_beginning, const int thread_zone_end) {
	if (verbose_flag) { Mutex.lock(); }
	for (int i = thread_zone_beginning; i < thread_zone_end; i++) {
		if (verbose_flag) {
			cout << LINE << " Worker " << thread_id << " START " << LINE;
		}
		
		BCH_objects[i] = make_shared<BCH_code_short_t2>(bitset<n>(vector_of_message_polynomials[i].to_string()));

		BCH_objects[i]->encode_bch();

		BCH_objects[i]->introduce_errors();

		BCH_objects[i]->print_original_codeword_and_received_codeword();

		BCH_objects[i]->Received_Data = bitset <k>(BCH_objects[i]->Received_Codeword.to_string().substr(0, k));

		bool status = BCH_objects[i]->decode_bch();

		// check the decoding status
		if (status == SUCCESS) {
			BCH_objects[i]->print_original_message_and_decoded_message();
			g_success_count++;
		} else {
			g_failure_count++;
		}

		if (verbose_flag) {
			cout << LINE << " Worker " << thread_id << " STOP *" << LINE << endl << endl;
		}
	}
	if (verbose_flag) { Mutex.unlock(); }
}

void populate_vector_of_message_polynomials (const unsigned char* str, const int thread_zone_beginning, const int thread_zone_end, 
									const int vector_zone_beginning, const int Bit_pos) {
	bitset<k> temp_bitset;
	int it = vector_zone_beginning;
	int bit_pos = Bit_pos;
	for (int i = thread_zone_beginning; i < thread_zone_end; i++) {
		for (int j = 0; j < 8; ++j) {
			temp_bitset[bit_pos] = str[i] & (1 << j);
			bit_pos++;
			if (bit_pos == k) {
				vector_of_message_polynomials[it] ^= temp_bitset;
				it++;
				temp_bitset.reset();
				bit_pos = 0;
			}
		}
	}
	if (bit_pos != 0) {
		vector_of_message_polynomials[it] ^= temp_bitset;
	}
}

void populate_unsigned_char_vectors (const int thread_zone_beginning, const int thread_zone_end, 
										const int vector_zone_beginning, const int Bit_pos) {
	int it = vector_zone_beginning;
	int bit_pos = Bit_pos;
	for (int i = thread_zone_beginning; i < thread_zone_end; i++) {
		for (int j = 0; j < 8; ++j) {
			recovered_charstream[i] ^= BCH_objects[it]->Decoded_Data[bit_pos] << j;
			modified_charstream[i] ^= BCH_objects[it]->Received_Data[bit_pos] << j;
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
		cout << "Invalid number of arguments"<< endl;
		return FAIL;
	}

	while ((c = getopt (argc, argv, "i:p:vh")) != FAIL) {
		switch (c) {
			case 'h':
				print_help_message(argv[0]);
				exit(0);
			case 'i':
				filename = optarg;
				fd = open(optarg, O_RDONLY);
				if (fd == FAIL) {
					cout << "Incorrect Image file name" << endl;
					return FAIL;
				}
				break;
			case 'p':
				error_probability = atoi(optarg);
				if (error_probability < 10 || error_probability > 10000000) {
					cout << "Error probability should be between (10 and 10000000)" << endl;
					return FAIL;
				}
				break;
			case 'v':
				verbose_flag = true;
				break;
			case '?':
					cout << optopt << endl;
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
		cout << "Argument parsing failed, type: \"" << argv[0] << " -h\" for more information" << endl;
		return 1;
	}
	
	struct stat sb;

	if (fstat(fd, &sb) == FAIL) {
		cout << "Failed to get file size\n";
		return 1;
	}

	unsigned char* buffer = (unsigned char*) mmap(
		NULL, 
		sb.st_size, 
		PROT_READ, MAP_PRIVATE, 
		fd, 
		0);

	string image_with_errors_path = filename.substr(0, filename.length()-4).append("_with_errors_BCH_code_short_t2.bmp");
	string image_fixed_path = filename.substr(0, filename.length()-4).append("_fixed_BCH_code_short_t2.bmp");

	remove(image_with_errors_path.c_str());
	remove(image_fixed_path.c_str());

	ofstream image_with_errors (image_with_errors_path.c_str(), ios::out | ios::app | ios::binary);
	ofstream image_fixed (image_fixed_path.c_str(), ios::out | ios::app | ios::binary);

	const ssize_t NUMBER_OF_MESSAGE_POLYNOMIALS = (sb.st_size*8)/k + 1;

	vector_of_message_polynomials.resize(NUMBER_OF_MESSAGE_POLYNOMIALS);
	BCH_objects.resize(NUMBER_OF_MESSAGE_POLYNOMIALS);
	recovered_charstream.resize(sb.st_size);
	modified_charstream.resize(sb.st_size);

	const int NUM_THREADS = thread::hardware_concurrency();
	cout << "number of detected threads: " << NUM_THREADS << endl;

	vector<thread> threads;

	const int THREAD_ZONE = NUMBER_OF_MESSAGE_POLYNOMIALS / NUM_THREADS;
	const int BIG_THREAD_ZONE = sb.st_size / NUM_THREADS;

	cout << "Parsing image file..." << endl;
	auto start = chrono::high_resolution_clock::now();

	int old_bit_pos = 0;
	int overhang = 0;

	for (int thread_id = 0; thread_id < NUM_THREADS; thread_id++) {
		int thread_zone_beginning = BIG_THREAD_ZONE*thread_id;
		int thread_zone_end = BIG_THREAD_ZONE*(thread_id + 1);

		if (thread_id == NUM_THREADS-1) {
			thread_zone_end = sb.st_size;
		}

		int bit_pos = (thread_zone_beginning*8) % k;

		if (bit_pos < old_bit_pos) { overhang++; }

		int vector_zone_beginning = THREAD_ZONE*thread_id + overhang;

		old_bit_pos = bit_pos;

		threads.push_back(thread(populate_vector_of_message_polynomials, buffer, thread_zone_beginning, 
									thread_zone_end, vector_zone_beginning, bit_pos));
	}

	for (auto& thread : threads) {
		thread.join();
	}

	auto stop = chrono::high_resolution_clock::now();
	auto dur = stop - start;
	auto duration = chrono::duration_cast<chrono::duration<float>>(dur);

	cout << "Parsing image file done and it took: ";
	cout << to_string(duration.count()).substr(0, to_string(duration.count()).length()-3) << " seconds" << endl;

	start = chrono::high_resolution_clock::now();

	BCH::initialize_BCH();

	threads.clear();

	if (!verbose_flag){
		cout<<"Please be patient, starting coding and decoding process..." << endl;
	}

	for (int thread_id = 0; thread_id < NUM_THREADS; thread_id++) {
		int thread_zone_beginning = THREAD_ZONE*thread_id;
		int thread_zone_end = THREAD_ZONE*(thread_id + 1);

		if (thread_id == NUM_THREADS-1) {
			thread_zone_end = NUMBER_OF_MESSAGE_POLYNOMIALS;
		}

		threads.push_back(thread(begin_main_process, thread_id, thread_zone_beginning, thread_zone_end));

	}

	for (auto& thread : threads) {
		thread.join();
	}

	stop = chrono::high_resolution_clock::now();
	dur = stop - start;
	auto main_duration = chrono::duration_cast<chrono::duration<float>>(dur);

	cout << "Coding and decoding process done and it took: ";
	cout << to_string(main_duration.count()).substr(0, to_string(main_duration.count()).length()-3) << " seconds" << endl;

	cout << "Converting modified and recovered data from bitsets to bytes..." << endl;

	start = chrono::high_resolution_clock::now();

	threads.clear();

	old_bit_pos = 0;
	overhang = 0;

	for (int thread_id = 0; thread_id < NUM_THREADS; thread_id++) {

		int thread_zone_beginning = BIG_THREAD_ZONE*thread_id;
		int thread_zone_end = BIG_THREAD_ZONE*(thread_id + 1);

		if (thread_id == NUM_THREADS-1) {
			thread_zone_end = sb.st_size;
		}

		int bit_pos = (thread_zone_beginning*8) % k;

		if (bit_pos < old_bit_pos) { overhang++; }

		int vector_zone_beginning = THREAD_ZONE*thread_id + overhang;

		old_bit_pos = bit_pos;

		threads.push_back(thread(populate_unsigned_char_vectors, thread_zone_beginning, thread_zone_end, 
									vector_zone_beginning, bit_pos));
	}

	for (auto& thread : threads) {
		thread.join();
	}

	stop = chrono::high_resolution_clock::now();
	dur = stop - start;
	duration = chrono::duration_cast<chrono::duration<float>>(dur);

	cout << "Converting modified and recovered data from bitsets to bytes done and it took: ";
	cout << to_string(duration.count()).substr(0, to_string(duration.count()).length()-3) << " seconds" << endl;

	for (int i = 0; i < HEADER_BYTES; i++) {
		image_with_errors << buffer[i];
		image_fixed << buffer[i];
	}

	for (int i = HEADER_BYTES; i < sb.st_size; i++) {
		image_with_errors << modified_charstream[i];
		image_fixed << recovered_charstream[i];
	}

	int difference_count = 0;
	for (int i = 0; i < NUMBER_OF_MESSAGE_POLYNOMIALS; i++) {
		difference_count += (vector_of_message_polynomials[i] ^ BCH_objects[i]->Decoded_Data).count();
	}

	// int total_fresh_diff = accumulate (BCH_objects.begin(), BCH_objects.end(), 0, 
	//  	[](int sum, shared_ptr<BCH_code_short_t2> const& obj) 
	// 	{ return sum + (obj->Codeword ^ obj->Received_Codeword).count(); });
	// if (total_fresh_diff != g_introduced_errors_count) { cout<< "Some of introduced errors cancelled each other out" << endl; }

	auto table_row_separator = []()->void { cout << setw(42) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl; };
	auto table_items = [](const string &lhs, const string &rhs)->void { cout << "| " << setw(40) << lhs << setw(20) << "| " + rhs << "|" << endl; };

	cout << endl << left;
	table_row_separator();
	table_items("Code used ", "(" + to_string(n) + "," + to_string(k) + "," + to_string(t*2+1 ) + ")");
	table_row_separator();
	table_items("Original image used ", filename);
	table_row_separator();
	table_items("Given error probability ", to_string(1/(float)error_probability));
	table_row_separator();
	table_items("Real error probability ", to_string(1/(float)((NUMBER_OF_MESSAGE_POLYNOMIALS*n)/g_introduced_errors_count)));
	table_row_separator();
	table_items("Number of all data bits ", to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * k));
	table_row_separator();
	table_items("Number of all data + redundant bits ", to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * n));
	table_row_separator();
	table_items("Number of all introduced errors ", to_string(g_introduced_errors_count));
	table_row_separator();
	table_items("Number of all message polynomials ", to_string(NUMBER_OF_MESSAGE_POLYNOMIALS));
	table_row_separator();
	table_items("Number of successful decodings", to_string(g_success_count));
	table_row_separator();
	table_items("Number of decoding errors ", to_string(g_failure_count));
	table_row_separator();
	table_items("Number of uncaught decoding errors ", to_string(g_uncaught_errors_count));
	table_row_separator();
	table_items("Number of over t errors in codewords ", to_string(g_big_errors_count));
	table_row_separator();
	table_items("Final data bits difference ", to_string(difference_count));
	table_row_separator();
	table_items("Encoding and decoding time ", to_string(main_duration.count()).substr(0, to_string(main_duration.count()).length()-3) + " seconds ");
	table_row_separator();

	cout << endl << "Results have been written to BCH_logs.txt file" << endl;
	cout << endl << "To view made images on linux use \"feh -F -Z --force-aliasing -d " << filename << " " << image_with_errors_path.c_str() << " " << image_fixed_path.c_str() << "\"" << endl << endl;

	ofstream BCH_logs;
	BCH_logs.open ("BCH_logs.txt", ios::out | ios::app);
	auto BCH_logs_items = [&](const string &lhs, const string &rhs)->void { BCH_logs << setw(37) << lhs << setw(20) << "| " + rhs << endl; };

	BCH_logs << endl << left;
	BCH_logs_items("Code used ", "(" + to_string(n) + "," + to_string(k) + "," + to_string(t*2+1 ) + ")");
	BCH_logs_items("Original image used ", filename);
	BCH_logs_items("Given error probability ", to_string(1/(float)error_probability));
	BCH_logs_items("Real error probability ", to_string(1/(float)((NUMBER_OF_MESSAGE_POLYNOMIALS*n)/g_introduced_errors_count)));
	BCH_logs_items("Number of all data bits ", to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * k));
	BCH_logs_items("Number of all data + redundant bits ", to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * n));
	BCH_logs_items("Number of all generated errors ", to_string(g_introduced_errors_count));
	BCH_logs_items("Number of all message polynomials ", to_string(NUMBER_OF_MESSAGE_POLYNOMIALS));
	BCH_logs_items("Number of successful decodings ", to_string(g_success_count));
	BCH_logs_items("Number of decoding errors ", to_string(g_failure_count));
	BCH_logs_items("Number of uncaught decoding errors", to_string(g_uncaught_errors_count));
	BCH_logs_items("Number of over t errors in codewords", to_string(g_big_errors_count));
	BCH_logs_items("Final data bits difference ", to_string(difference_count));
	BCH_logs_items("Encoding and decoding time ", to_string(main_duration.count()).substr(0, to_string(main_duration.count()).length()-3) + " seconds ");

	munmap(buffer, sb.st_size);
	close(fd);

	image_with_errors.close();
	image_fixed.close();
	BCH_logs.close();

	return 0;
}
