#include "bch4830.hpp"

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
			for (int root = 1; root < 2*t; root++)
				if (element == root) {
					root_found = true;
					break;
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

bitset<n> BCH_code_short_t3::encode_bch(const bitset<n> &Data) {
	bitset<n> Codeword;
	bitset<n> Shifted_Data = Data<<(n-k);
	bitset<n> rb = divide_bitset_polynomials(Shifted_Data, generator_polynomial_bitset).first;
	Codeword = Shifted_Data ^ rb;
	if (verbose) {
		cout << endl;
	}
	return Codeword;
}

vector <int> BCH_code_short_t3::calculate_syndromes(const bitset<n> &Received_Codeword, bool &errors_in_codeword) {
	vector <int> syndromes(2*t+1);
	for (int i = 1; i <= 2*t; i++) {
		syndromes[i] = 0;
		for (int j = 0; j < n; j++) {
			if (Received_Codeword[n - j - 1] != 0) {
				syndromes[i] ^= alpha_to[(i*j) % GFB];
			}
		}
		syndromes[i] = index_of[syndromes[i]];
		if (syndromes[i] != -1) {
			errors_in_codeword = true;
		}
	}
	if (verbose) {
		cout << "syndromes= ( ";
		for (int i = 1; i <= 2*t; i++) {
			cout << syndromes[i] << " ";
		}
		cout << ")" << endl;
	}
	return syndromes;
}

pair<bitset<k>, bool> BCH_code_short_t3::decode_bch(const bitset <n> &Received_Codeword) {
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
	vector <int> syndromes = calculate_syndromes(Received_Codeword, errors_in_codeword);
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
			if (verbose) {
				cout << "Sigma(x) = ";
			}
			for (int i = 0; i <= l[u]; i++) {
				elp[u][i] = index_of[elp[u][i]];
				if (verbose) {
					cout << elp[u][i] << " ";
				}
			}
			// Find roots of the error location polynomial:
			// Chien search
			if (verbose) {
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
					if (verbose) {
						cout << error_locations.back() << " ";
					}
				}
			}
			if (verbose) {
				cout << endl;
			}
			if (error_locations.size() == (unsigned)l[u]) {
			// no. roots = degree of elp hence <= t errors
				for (auto const &error_location : error_locations) {
					auto err_loc = (n-1-error_location + n) % n;
					if (err_loc < 0 || err_loc >= n) {
						if (verbose) {
							cout<<endl<<"Incomplete decoding: errors detected (err_loc out of bound: err_loc = " << err_loc << ")"<<endl;
						}
						bitset<k> Decoded_Message(Decoded_Codeword.to_string().substr(0, k));
						return {Decoded_Message, false};
					} else {
						Decoded_Codeword.flip(err_loc);
					}
				}
			} else {// elp has degree > t hence cannot solve
				if (verbose) {
					cout<<endl<<"Incomplete decoding: errors detected (elp has degree > t) "<<endl;
				}
				bitset<k> Decoded_Message(Decoded_Codeword.to_string().substr(0, k));
				return {Decoded_Message, false};
			}
		} else {
			if (verbose) {
				cout<<endl<<"Incomplete decoding: errors detected (l[u] > t: l[u] = " << l[u] << ")"<<endl;
			}
			bitset<k> Decoded_Message(Decoded_Codeword.to_string().substr(0, k));
			return {Decoded_Message, false};
		}
	} else {
		if (verbose) {
			cout << "No errors found" << endl;
		}
	}
	bitset<k> Decoded_Message(Decoded_Codeword.to_string().substr(0, k));
	return {Decoded_Message, true};
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

pair<bitset <n>, bitset <n>> BCH_code_short_t3::divide_bitset_polynomials(const bitset <n> &Dividend, const bitset <n> &Divisor) {
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

void BCH_code_short_t3::print_original_codeword_and_received_codeword(const bitset <n> &Codeword, const bitset <n> &Received_Codeword) {
	if (verbose) {
		cout << "c(x) = "<< Codeword << endl;
		cout << "r(x) = "<< Received_Codeword << endl;
	}
	bitset <n> bit_difference = Codeword ^ Received_Codeword;
	if (bit_difference.count()) {
		if (verbose) {
			cout<<"       "<<(bit_difference.to_string(' ','^'))<<endl;
			if (bit_difference.count() == 1) {
				cout << "Position of " << bit_difference.count() << " error in the received codeword at ^." << endl;
			} else {
				cout << "Positions of " << bit_difference.count() << " errors in the received codeword at ^." << endl;
			}
		}
		total_errors += bit_difference.count();
		if (bit_difference.count() > t) {
			big_errors++;
		}
	} else {
		if (verbose) {
			cout << "No errors in the received codeword." << endl;
		}
	}
}

void BCH_code_short_t3::print_original_message_and_decoded_message(const bitset <n> &Data, const bitset <k> &Decoded_Data) {
	if (verbose) {
		cout <<  DASH_LINE << "Results" << DASH_LINE << endl;
	}
	bitset<k> Original_Data(Data.to_string().substr(n-k, n));
	if (verbose) {
		cout << "Original Data  = " << Original_Data << endl;
		cout << "Recovered Data = " << Decoded_Data<< endl;
	}
	bitset <k> bit_difference = Original_Data ^ Decoded_Data;
	if (bit_difference.count()) {
		if (verbose) {
			cout << "                 " << (bit_difference.to_string(' ','^')) << endl;
			if (bit_difference.count() == 1) {
				cout << "Position of " << bit_difference.count() << " message decoding error at ^." << endl;
			} else {
				cout << "Positions of " << bit_difference.count() << " message decoding errors at ^." << endl;
			}
		}
		uncaught_errors++;
	} else {
		if (verbose) {
			cout <<"Succesful decoding." << endl;
		}
	}
}

bitset <n> BCH_code_short_t3::introduce_errors(const bitset <n> &Codeword, const int error_probability) {
	bitset <n> Received_Codeword = Codeword;
	random_device rd;
	mt19937 gen (rd());
	uniform_int_distribution<> d (0, error_probability);
	for (int i = n - 1; i >= 0; i--) {
		int randnum = d(gen);
		if (randnum == 0) {
			Received_Codeword.flip(i);
		}
	}
	return Received_Codeword;
}

void print_help_message(const char *file_name) {
	cout << "Usage:\n"
		<< file_name << " [-h] -i <image> -p err_prob [-v]\n\n"
		<< "Options:\n"
		<< "  -i <file_name>, --image <file_name>	Choose one image from images folder, example: images/image2.bmp.\n"
		<< "  -p err_prob, --probability err_prob	Give probability between (1 and 10000000) that an error will occur in \n"
		<< "					  the codeword during a simulated transmission through a noisy medium, example: 1000.\n"
		<< "Optional arguments:\n"
		<< "  -h, --help	Show this help message.\n"
		<< "  -v, --verbose	Enable verbose encoding and decoding logs which will print out the whole process to the \n"
		<< "		  terminal, is disabled by default. WARNING! This option causes the threads to run sequenitally instead \n"
		<< " 		  of in parallel which combined with printing operations to console causes a severe performance degradation.\n";
}

void begin_main_process(const int thread_id, const int thread_zone_beginning, const int thread_zone_end) {
	if (verbose) { Mutex.lock(); }
	for (int i = thread_zone_beginning; i < thread_zone_end; i++) {
		if (verbose) {
			cout << LINE << " Worker " << thread_id << " START " << LINE;
		}
		BCH_objects[i] = make_unique<BCH_code_short_t3>(bitset<n>(vector_of_bitsets[i].to_string()), error_probability);
		vector_of_modified_bitsets[i] = BCH_objects[i] -> Received_Data;
		vector_of_recovered_bitsets[i] = BCH_objects[i] -> Decoded_Data.first;
		if (verbose) {
			cout << LINE << " Worker " << thread_id << " STOP *" << LINE << endl << endl;
		}
	}
	if (verbose) { Mutex.unlock(); }
}

void populate_vector_of_bitsets (const unsigned char* str, const int thread_zone_beginning, const int thread_zone_end, 
									const int vector_zone_beginning, const int Bit_pos) {
	bitset<k> temp_bitset;
	int it = vector_zone_beginning;
	int bit_pos = Bit_pos;
	for (int i = thread_zone_beginning; i < thread_zone_end; i++) {
		for (int j = 0; j < 8; ++j) {
			temp_bitset[bit_pos] = str[i] & (1 << j);
			bit_pos++;
			if (bit_pos == k) {
				vector_of_bitsets[it] ^= temp_bitset;
				it++;
				temp_bitset.reset();
				bit_pos = 0;
			}
		}
	}
	if (bit_pos != 0) {
		vector_of_bitsets[it] ^= temp_bitset;
	}
}

void populate_unsigned_char_vectors (const int thread_zone_beginning, const int thread_zone_end, 
										const int vector_zone_beginning, const int Bit_pos) {
	int it = vector_zone_beginning;
	int bit_pos = Bit_pos;
	for (int i = thread_zone_beginning; i < thread_zone_end; i++) {
		for (int j = 0; j < 8; ++j) {
			recovered_charstream[i] ^= vector_of_recovered_bitsets[it][bit_pos] << j;
			modified_charstream[i] ^= vector_of_modified_bitsets[it][bit_pos] << j;
			bit_pos++;
			if (bit_pos == k) {
				it++;	
				bit_pos = 0;
			}
		}
	}
}

int main(int argc, const char* argv[]) {
	int fd;
	if (argc != 2 && argc != 5 && argc != 6) {
		cout << "Invalid number of arguments, type: \"" << argv[0] << " -h\" for more information" << endl;
		return 0;
	}

	if (strcasecmp(argv[1], "-h") == 0 || strcasecmp(argv[1], "--help") == 0) {
		print_help_message(argv[0]);
		return 0;
	} else if (strcasecmp(argv[1], "-i") == 0 || strcasecmp(argv[1], "--image") == 0) {
		fd = open(argv[2], O_RDONLY);
		if (fd == -1) {
			cout << "Incorrect Image file name, type: \"" << argv[0] << " -h\" for more information" << endl;
			return 0;
		}
	} else {
		cout << "Invalid argument, type: \"" << argv[0] << " -h\" for more information" << endl;
		return 0;
	}

	if (strcasecmp(argv[3], "-p") == 0 || strcasecmp(argv[3], "--probability") == 0) {
		error_probability = stoi(argv[4]);
		if (error_probability < 1 || error_probability > 10000000) {
			cout << "Error probability should be between (1 and 10000000), type: \""
			<< argv[0] << " -h\" for more information" << endl;
			return 0;
		}
	} else {
		cout << "Invalid argument, type: \"" << argv[0] << " -h\" for more information" << endl;
		return 0;
	}

	if (argc == 6) {
		if (strcasecmp(argv[5], "-v") == 0 || strcasecmp(argv[5], "--verbose") == 0) {
			verbose = true;
		} else {
			cout << "Invalid argument, type: \"" << argv[0] << " -h\" for more information" << endl;
			return 0;
		}
	}

	struct stat sb;
	if (fstat(fd, &sb) == -1) {
		cout << "Failed to get file size\n";
		return 1;
	}
	unsigned char* buffer = (unsigned char*) mmap(
		NULL, 
		sb.st_size, 
		PROT_READ, MAP_PRIVATE, 
		fd, 
		0);

	string image_with_errors_path = string(argv[2]).substr(0, string(argv[2]).length()-4).append("_with_errors_BCH_code_short_t3.bmp");
	string image_fixed_path = string(argv[2]).substr(0, string(argv[2]).length()-4).append("_fixed_BCH_code_short_t3.bmp");
	remove(image_with_errors_path.c_str());
	remove(image_fixed_path.c_str());
	ofstream image_with_errors (image_with_errors_path.c_str(), ios::out | ios::app | ios::binary);
	ofstream image_fixed (image_fixed_path.c_str(), ios::out | ios::app | ios::binary);
	const ssize_t NUMBER_OF_MESSAGE_POLYNOMIALS = (sb.st_size*8)/k + 1;
	vector_of_bitsets.resize(NUMBER_OF_MESSAGE_POLYNOMIALS);
	vector_of_modified_bitsets.resize(NUMBER_OF_MESSAGE_POLYNOMIALS);
	vector_of_recovered_bitsets.resize(NUMBER_OF_MESSAGE_POLYNOMIALS);
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
		threads.push_back(thread(populate_vector_of_bitsets, buffer, thread_zone_beginning, 
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
	if (!verbose){
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
		difference_count += (vector_of_bitsets[i] ^ vector_of_recovered_bitsets[i]).count();
	}
	cout << endl << left << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Code used " << setw(20) << "| (" + to_string(n) + "," + to_string(k) + "," + to_string(t*2+1 ) + ")" << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Original image used " << setw(20) << "| " + string(argv[2]) << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Given error probability " << setw(20) << "| 1 in " + to_string(error_probability) << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Real error probability " << setw(20) << "| 1 in " + to_string((NUMBER_OF_MESSAGE_POLYNOMIALS*n)/total_errors) << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of all bits " << setw(20) << "| " + to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * k) << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of all bits + redundant bits " << setw(20) << "| " + to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * n) << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of all generated errors " << setw(20) << "| " + to_string(total_errors) << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of all message polynomials " << setw(20) << "| " + to_string(NUMBER_OF_MESSAGE_POLYNOMIALS) << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of successful decodings " << setw(20) << "| " + to_string(success) << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of decoding errors " << setw(20) << "| " + to_string(failure) << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of uncaught decoding errors " << setw(20) << "| " + to_string(uncaught_errors) << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Number of over t errors in codeword " << setw(20) << "| " + to_string(big_errors) << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Final bit difference " << setw(20) << "| " + to_string(difference_count) << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	cout << "|" << setw(37) <<" Encoding and decoding time " << setw(20) << "| " + to_string(main_duration.count()).substr(0, to_string(main_duration.count()).length()-3) + " seconds " << "|" << endl;
	cout << setw(38) << setfill('-') << "+" << setw(20) << "+" << setw(1) << "+" << setfill (' ') << endl;
	
	cout << endl << "Results have been written to BCH_logs.txt file" << endl;
	cout << endl << "To view made images on linux use \"feh -F -Z --force-aliasing -d " << argv[2] << " " << image_with_errors_path.c_str() << " " << image_fixed_path.c_str() << "\"" << endl << endl;

	ofstream BCH_logs;
	BCH_logs.open ("BCH_logs.txt", ios::out | ios::app);
	BCH_logs << left << endl << setw(37) <<" Code used " << setw(20) << "| (" + to_string(n) + "," + to_string(k) + "," + to_string(t*2+1 ) + ")" << endl;
	BCH_logs << setw(37) <<" Original image used " << setw(20) << "| " + string(argv[2]) << endl;
	BCH_logs << setw(37) <<" Given error probability " << setw(20) << "| 1 in " + to_string(error_probability) << endl;
	BCH_logs << setw(37) <<" Real error probability " << setw(20) << "| 1 in " + to_string((NUMBER_OF_MESSAGE_POLYNOMIALS*n)/total_errors) << endl;
	BCH_logs << setw(37) <<" Number of all bits " << setw(20) << "| " + to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * k) << endl;
	BCH_logs << setw(37) <<" Number of all bits + redundant bits " << setw(20) << "| " + to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * n) << endl;
	BCH_logs << setw(37) <<" Number of all generated errors " << setw(20) << "| " + to_string(total_errors) << endl;
	BCH_logs << setw(37) <<" Number of all message polynomials " << setw(20) << "| " + to_string(NUMBER_OF_MESSAGE_POLYNOMIALS) << endl;
	BCH_logs << setw(37) <<" Number of successful decodings " << setw(20) << "| " + to_string(success) << endl;
	BCH_logs << setw(37) <<" Number of decoding errors " << setw(20) << "| " + to_string(failure) << endl;
	BCH_logs << setw(37) <<" Number of uncaught decoding errors " << setw(20) << "| " + to_string(uncaught_errors) << endl;
	BCH_logs << setw(37) <<" Number of over t errors in codeword " << setw(20) << "| " + to_string(big_errors) << endl;
	BCH_logs << setw(37) <<" Final bit difference " << setw(20) << "| " + to_string(difference_count) << endl;
	BCH_logs << setw(37) <<" Encoding and decoding time " << setw(20) << "| " + to_string(main_duration.count()).substr(0, to_string(main_duration.count()).length()-3) + " seconds " << endl;

	munmap(buffer, sb.st_size);
	close(fd);
	image_with_errors.close();
	image_fixed.close();
	BCH_logs.close();
	return 0;
}
