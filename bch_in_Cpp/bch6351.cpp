#include "bch6351.hpp"
#include "bch_logger.hpp"

// use overload pattern here
// template<typename... Ts> struct overload : public Ts... { using Ts::operator()...; };
// template<class... Ts> overload(Ts...) -> overload<Ts...>;

using bchType = std::variant<std::unique_ptr<Bch6351>, std::unique_ptr<Bch6345>, 
                             std::unique_ptr<Bch4836>, std::unique_ptr<Bch4830>>;


template<typename T>
// (CURRENT VARIANT BCH CLASS) Alias for getting the actual type of the Bch code type class used in variant
using CVBC = std::remove_pointer_t<std::decay_t<T>>;

std::vector<std::bitset <Bch6351::k_>> Bch6351::vector_of_message_polynomials;
std::vector<std::bitset <Bch6345::k_>> Bch6345::vector_of_message_polynomials;
std::vector<std::bitset <Bch4836::k_>> Bch4836::vector_of_message_polynomials;
std::vector<std::bitset <Bch4830::k_>> Bch4830::vector_of_message_polynomials;

static BchLogger bch_logger;
static globalCounters counters{};

template <size_t X>
int bch::MSB(
		const std::bitset <X> &Polynomial)
{
	return (X + (GFB-X) - std::countl_zero(Polynomial.to_ullong()));
}

template <size_t N>
void bch::verbosePolynomial(
		const std::bitset <N> &Polynomial)
{
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

/**
 * Read the primitive polynomial of degree 6 from binary form
*/
template <class bch_class>
void readPrimitivePolynomial(
		mathHelper<bch_class::n_>* bch_math)
{
	bch_math->primitive_polynomial_int = bch_math->primitive_polynomial_bitset.to_ulong();
	std::cout << "Primitive polynomial:" << std::endl << "primitive_polynomial_bitset(x) = ";
	bch::verbosePolynomial(bch_math->primitive_polynomial_bitset);
}

// done
template<class bch_class>
void generateGaloisField(
		mathHelper<bch_class::n_>* bch_math)
{
	bch_math->index_of[0] = -1;
	bch_math->m = bch::MSB(bch_math->primitive_polynomial_bitset);
	for (int i = 0; i < GFB; i++) {
		if (i < bch_math->m) {
			bch_math->alpha_to[i] = 1 << i;
			bch_math->index_of[bch_math->alpha_to[i]] = i;
		} else {
			if (bch_math->alpha_to[i - 1] >= 32) {
				bch_math->alpha_to[i] = (bch_math->alpha_to[i - 1] << 1) ^ 
										bch_math->primitive_polynomial_int;
			} else {
				bch_math->alpha_to[i] = bch_math->alpha_to[i - 1] << 1;
			}
			bch_math->index_of[bch_math->alpha_to[i]] = i;
		}
	}
}

template <class bch_class>
void generateGeneratorPolynomial(
		mathHelper<bch_class::n_>* bch_math)
{
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
		for (int i = 1; i < bch_math->m; i++) {
			coset.push_back((coset[i-1] << 1) % GFB);
			// TODO: maybe implement if initializer
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
			for (int root = 1; root <= 2*bch_class::t_; root++) {
				if (element == root) {
					root_found = true;
					break;
				}
			}
			if(root_found) {
				bch_math->zeros_cosets.push_back(coset);
				break; 
			}
		}
		if (bch_math->zeros_cosets.size() == bch_class::t_) {
			break;
		}
	}
	std::vector <int> min_polynomials;
	for (const auto& zero_coset : bch_math->zeros_cosets) {
		int product = bch_math->alpha_to[zero_coset[0]] ^ 2; // (ax + x)
		for (uint i = 1; i < zero_coset.size(); i++) {
			product = bch::multiplyIntPolynomials(product, bch_math->alpha_to[zero_coset[i]] ^ 2);
		}
		product %= GFB+1;
		product ^= bch_math->primitive_polynomial_int;
		min_polynomials.push_back(product);
	}
	int generator_polynomial = min_polynomials[0];
	for (uint i=1; i<bch_class::t_; i++) {
		generator_polynomial = bch::multiplyIntPolynomials(generator_polynomial, min_polynomials[i]);
	}
	bch_math->generator_polynomial_bitset = generator_polynomial;
	bch::reverseBitset(bch_math->generator_polynomial_bitset, bch_class::k_-1);
	std::cout << "This is a (" << bch_class::n_ << "," << bch_class::k_ << "," << bch_class::t_*2+1 << ") binary bch code" << std::endl;
	std::cout << "g(x) is " << bch_math->generator_polynomial_bitset.to_string().substr(bch_class::n_ - (bch_class::n_-bch_class::k_+1)) << std::endl;
}

template<class bch_class>
void encodeBch(
		polynomialData<bch_class::n_, bch_class::k_>& polynomials, 
		const mathHelper<bch_class::n_>* bch_math)
{	
	polynomials.encoded.codeword <<= (bch_class::n_-bch_class::k_);
	std::bitset<bch_class::n_> rb = bch::divideBitsetPolynomials(polynomials.encoded.codeword, bch_math->generator_polynomial_bitset).first;
	polynomials.encoded.codeword ^= rb;
	bch::reverseBitset(polynomials.encoded.codeword, 0);
	bch_logger.log("\n");
}

template<class bch_class>
std::vector <int> calculateSyndromes(
		const polynomialData<bch_class::n_, bch_class::k_>& polynomials, 
		const mathHelper<bch_class::n_>* bch_math, 
		bool &errors_in_codeword)
{
	std::vector <int> syndromes(2*bch_class::t_+1);
	for (int i = 1; i <= 2*bch_class::t_; i++) {
		syndromes[i] = 0;
		for (int j = 0; j < bch_class::n_; j++) {
			if (polynomials.received.codeword[j] != 0) {
				syndromes[i] ^= bch_math->alpha_to[(i*j) % GFB];
			}
		}
		syndromes[i] = bch_math->index_of[syndromes[i]];
		// when there are no errors all syndromes should be -1
		if (syndromes[i] != -1) {
			errors_in_codeword = true;
		}
	}

	bch_logger.log("syndromes= ( ");
		for (int i = 1; i <= 2*bch_class::t_; i++) {
			bch_logger.log(syndromes[i], " ");
		}
	bch_logger.log(")\n");

	return syndromes;
}

template<class bch_class>
status decodeBch(
		polynomialData<bch_class::n_, bch_class::k_>& polynomials,
		const mathHelper<bch_class::n_>* bch_math) 
{
	/*
	 Simon Rockliff's implementation of Berlekamp's algorithm.
	 Assume we have received bits in Received_Codeword[i], i=0..(bch_class::n_-1).
	
	 Compute the 2*bch_class::t_ syndromes by substituting alpha^i into rec(X) and
	 evaluating, storing the syndromes in syndromes[i], i=1..2t (leave syndromes[0] zero) .
	 Then we use the Berlekamp algorithm to find the error location polynomial
	 elp[i].
	
	 If the degree of the elp is >bch_class::t_, then we cannot correct all the errors, and
	 we have detected an uncorrectable error pattern. We output the information
	 bits uncorrected.
	
	 If the degree of elp is <=bch_class::t_, we substitute alpha^i , i=1..bch_class::n_ into the elp
	 to get the roots, hence the inverse roots, the error location numbers.
	 This step is usually called "Chien's search".
	
	 If the number of errors located is not equal the degree of the elp, then
	 the decoder assumes that there are more than bch_class::t_ errors and cannot correct
	 them, only detect them. We output the information bits uncorrected.
	*/

	bool errors_in_codeword = false;
	std::vector <int> syndromes = calculateSyndromes<bch_class>(polynomials, bch_math, errors_in_codeword);
	polynomials.decoded.codeword = polynomials.received.codeword;

	// get them back to normal:
	bch::reverseBitset(polynomials.encoded.codeword, 0);
	bch::reverseBitset(polynomials.received.codeword, 0);
	if (errors_in_codeword) {
		/*
		 Compute the error location polynomial via the Berlekamp
		 iterative algorithm. Following the terminology of Lin and
		 Costello's book :   d[u] is the 'mu'th discrepancy, where
		 u='mu'+1 and 'mu' is the step number
		 ranging from -1 to 2*bch_class::t_ (see L&C),  l[u] is the degree of
		 the elp at that step, and u_l[u] is the difference between
		 the step number and the degree of the elp. 
		*/

		// WARNING
		// THIS CODE MUMBO JUMBO CAN MAKE YOU NAUSEOUS
		// WARNING

		std::vector <int> error_locations;
		int elp[100][100] = {{0, 0}}; // error location polynomial
		int d[100] = {0}, l[100] = {0}, u_lu[100] = {0};
		for (int i = 1; i < bch_class::t_*2; i++) {
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
					elp[u][i] = bch_math->index_of[elp[u][i]];
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
				for (int i = 0; i < 2*bch_class::t_; i++) {
					elp[u + 1][i] = 0;
				}
				for (int i = 0; i <= l[q]; i++) {
					if (elp[q][i] != -1) {
						elp[u + 1][i + u - q] = bch_math->alpha_to[(d[u] + 
							GFB - d[q] + elp[q][i]) % GFB];
					}
				}
				for (int i = 0; i <= l[u]; i++) {
					elp[u + 1][i] ^= elp[u][i];
					elp[u][i] = bch_math->index_of[elp[u][i]];
				}
			}
			u_lu[u + 1] = u - l[u + 1];
			// form (u+1)th discrepancy
			if (u < 2*bch_class::t_) {
			// no discrepancy computed on last iteration
				if (syndromes[u + 1] != -1) {
					d[u + 1] = bch_math->alpha_to[syndromes[u + 1]];
				} else {
				d[u + 1] = 0;
				}
				for (int i = 1; i <= l[u + 1]; i++) {
					if ((syndromes[u + 1 - i] != -1) && (elp[u + 1][i] != 0)) {
						d[u + 1] ^= bch_math->alpha_to[(syndromes[u + 1 - i] + 
							bch_math->index_of[elp[u + 1][i]]) % GFB];
					}
				}
				// put d[u+1] into index form
				d[u + 1] = bch_math->index_of[d[u + 1]];
			}
		} while ((u < 2*bch_class::t_) && (l[u + 1] <= bch_class::t_));
		u++;
		if (l[u] <= bch_class::t_) {// Can correct errors
			// put elp into index form
			bch_logger.log("Sigma(x) = ");
			for (int i = 0; i <= l[u]; i++) {
				elp[u][i] = bch_math->index_of[elp[u][i]];
				bch_logger.log(elp[u][i], " ");
			}
			// Find roots of the error location polynomial:
			// Chien search
			bch_logger.log("\nRoots: ");
			for (int i = 1; i <= GFB; i++) {
				int q = 1;
				for (int j = 1; j <= l[u]; j++) {
					elp[u][j] = (elp[u][j] + j) % GFB;
					q ^= bch_math->alpha_to[elp[u][j]];
				}
				// Store error location number indices
				if (!q) {
					error_locations.push_back(i  - 1 -(GFB-bch_class::n_));
					bch_logger.log(error_locations.back(), " ");
				}
			}
			bch_logger.log("\n");
			if (error_locations.size() == static_cast<unsigned>(l[u])) {
			// no. roots = degree of elp hence <= bch_class::t_ errors
				for (auto const &error_location : error_locations) {
					auto err_loc = (bch_class::n_-1-error_location + bch_class::n_) % bch_class::n_;
					if (err_loc < 0 || err_loc >= bch_class::n_) {
						bch_logger.log("\nIncomplete decoding: errors detected (err_loc out of bound: err_loc = ", err_loc, ")\n");
						// we can skip this thanks to union
						// polynomials.decoded.data = (std::bitset<bch_class::k_>)Decoded_Codeword.to_string().substr(0, bch_class::k_);
						bch::reverseBitset(polynomials.decoded.codeword, 0);
						return FAIL;
					} else {
						polynomials.decoded.codeword.flip(err_loc);
					}
				}
			} else {// elp has degree > bch_class::t_ hence cannot solve
				bch_logger.log("\nIncomplete decoding: errors detected (elp has degree > bch_class::t_)\n");
				// we can skip this thanks to union
				// polynomials.decoded.data = (std::bitset<bch_class::k_>)Decoded_Codeword.to_string().substr(0, bch_class::k_);
				bch::reverseBitset(polynomials.decoded.codeword, 0);
				return FAIL;
			}
		} else {
			bch_logger.log("\nIncomplete decoding: errors detected (l[u] > bch_class::t_: l[u] = ", l[u], ")\n");
			// we can skip this thanks to union
			// polynomials.decoded.data = (std::bitset<bch_class::k_>)Decoded_Codeword.to_string().substr(0, bch_class::k_);
			bch::reverseBitset(polynomials.decoded.codeword, 0);
			return FAIL;
		}
	} else {
		bch_logger.log("No errors found\n");
	}
	// we can skip this thanks to union
	// polynomials.decoded.data = (std::bitset<bch_class::k_>)Decoded_Codeword.to_string().substr(0, bch_class::k_);
	bch::reverseBitset(polynomials.decoded.codeword, 0);
	return SUCCESS;
}

template<size_t X>
void bch::reverseBitset(
		std::bitset <X> &Polynomial, 
		int Shift)
{
    for(size_t i = 0; i < X/2; ++i) {
    	bool temp_bit = Polynomial[i];
    	Polynomial[i] = Polynomial[X-i-1];
    	Polynomial[X-i-1] = temp_bit;
    }
	Polynomial >>= Shift;
}

template<size_t X>
std::bitset <X> bch::tempReverseBitset(
		std::bitset <X> Polynomial, 
		int Shift)
{
    for(size_t i = 0; i < X/2; ++i) {
    	bool temp_bit = Polynomial[i];
    	Polynomial[i] = Polynomial[X-i-1];
    	Polynomial[X-i-1] = temp_bit;
    }
	return (Polynomial >>= Shift);
}

template <size_t X>
std::pair<std::bitset <X>, std::bitset <X>> bch::divideBitsetPolynomials(
		const std::bitset <X> &dividend, 
		const std::bitset <X> &divisor)
{
	std::bitset <X> quotient, remainder = dividend;
 	while (bch::MSB(remainder) >= bch::MSB(divisor)) {
		int shift = bch::MSB(remainder) - bch::MSB(divisor);
		remainder ^= divisor << shift;
		quotient.flip(shift); 
	}
	return {remainder, quotient};
}

int bch::multiplyIntPolynomials(
		int mulitplicand, 
		int multiplicator)
{
	int product = 0;
	while (mulitplicand > 0) {
		if (mulitplicand & 1) {
			product ^= multiplicator;
		}
		multiplicator <<= 1;
		mulitplicand >>= 1;
	}
	return product;
}

template<class bch_class>
void compareAndPrintCodewords(
		const polynomialData<bch_class::n_, bch_class::k_>& polynomials)
{	
	bch_logger.log("c(x) = ", bch::tempReverseBitset(polynomials.encoded.codeword, 0), "\n");
	bch_logger.log("r(x) = ", bch::tempReverseBitset(polynomials.received.codeword, 0), "\n");

	std::bitset <bch_class::n_> bit_difference = polynomials.encoded.codeword ^ polynomials.received.codeword;
	bch::reverseBitset(bit_difference, 0);

	if (bit_difference.count() != 0) {
		if (bch_logger.enable_logging_) {
			bch_logger.log("       ", (bit_difference.to_string(' ','^')), "\n");
			bch_logger.log("Position of ", bit_difference.count());
			bit_difference.count() == 1? bch_logger.log(" error "): bch_logger.log(" errors ");
			bch_logger.log("in the received codeword at ^.\n");
		}
		if (bit_difference.count() > static_cast<unsigned>(bch_class::t_)) {
			counters.big_errors_count++;
		}
	} else {
		bch_logger.log("No errors in the received codeword.\n");
	}
}

template<class bch_class>
void compareAndPrintData(
		const polynomialData<bch_class::n_, bch_class::k_>& polynomials)
{
	bch_logger.log(DASH_LINE, "Results", DASH_LINE,"\n");

	bch_logger.log("Original Data  = ", polynomials.encoded.data, "\n");
	bch_logger.log("Recovered Data = ", polynomials.decoded.data, "\n");

	std::bitset <bch_class::k_> bit_difference = polynomials.encoded.data ^ polynomials.decoded.data;

	if (bit_difference.count() != 0) {
		bch_logger.log("                 ", (bit_difference.to_string(' ','^')), "\n");
		bch_logger.log("Position of ", bit_difference.count()," message decoding");
		bit_difference.count() == 1? bch_logger.log(" error "): bch_logger.log(" errors ");
		bch_logger.log("at ^.\n");
		counters.uncaught_errors_count++;
	} else {
		bch_logger.log("Successful decoding.\n");
	}
}

template<class bch_class>
void introduceErrors(
		polynomialData<bch_class::n_, bch_class::k_>& polynomials) 
{
	polynomials.received.codeword = polynomials.encoded.codeword;
	std::random_device rd;
	std::mt19937 gen (rd());
	std::uniform_int_distribution<> d (0, bch::error_probability);
	for (int i = 0; i < bch_class::n_; i++) {
		int randnum = d(gen);
		if (randnum == 0) {
			counters.introduced_errors_count ++;
			polynomials.received.codeword.flip(i);
		}
	}
}

void printHelpMessage(
		const char *file_name)
{
	std::cout << "Usage:\n"
		<< file_name << " [-h] -i image -p err_prob -c code_type [-v]\n\n"
		<< "Options:\n"
		<< "  -i image			Choose one image from images folder, example: images/image2.bmp.\n"
		<< "  -p err_prob		Give probability between (10 and 10000000) that a 1 in err_prob error will occur in \n"
		<< "			   			the codeword during a simulated transmission through a noisy medium, example: 1000.\n"
		<< "  -c code_type		Choose bch code type:\n"
		<< "						0 - BCH6351\n"
		<< "						1 - BCH6345\n"
		<< "						2 - BCH4836\n"
		<< "						3 - BCH4830\n"
		<< "Optional arguments:\n"
		<< "  -h	Show this help message.\n"
		<< "  -v	Enable verbose_flag encoding and decoding logs which will print out the whole process to the \n"
		<< "	   terminal, is disabled by default. WARNING! This option causes the threads to run sequentially instead \n"
		<< " 	   of in parallel which combined with printing operations to console causes a severe performance degradation.\n";
}

void beginMainProcess(
		std::vector<bchType>& objs,
		const int thread_id, 
		const threadZones& zone)
{	
	if (bch_logger.enable_logging_) { bch::Mutex.lock(); }
	for (int i = zone.MPTG_beginning; i < zone.MPTG_end; i++) {
		bch_logger.log(LINE, " Worker ", thread_id," START ", LINE);

		switch (bch::code_type) {
			case (BCH6351):
				objs[i] = std::make_unique<Bch6351>(Bch6351::vector_of_message_polynomials[i]);
				break;
			case (BCH6345):
				objs[i] = std::make_unique<Bch6345>(Bch6345::vector_of_message_polynomials[i]);
				break;
			case (BCH4836):
				objs[i] = std::make_unique<Bch4836>(Bch4836::vector_of_message_polynomials[i]);
				break;
			case (BCH4830):
				objs[i] = std::make_unique<Bch4830>(Bch4830::vector_of_message_polynomials[i]);
				break;
		}

		std::visit([&](auto& obj) -> void {
				using current_class = CVBC<decltype(obj.get())>;

				encodeBch<current_class>(obj->codeword_polynomials_, bch::bch_math<current_class::n_>);

				introduceErrors<current_class>(obj->codeword_polynomials_);

				compareAndPrintCodewords<current_class>(obj->codeword_polynomials_);

				if (status st = decodeBch<current_class>(obj->codeword_polynomials_, bch::bch_math<current_class::n_>); st == SUCCESS) {
					compareAndPrintData<current_class>(obj->codeword_polynomials_);
					counters.success_count++;
				} else {
					counters.failure_count++;
				}
		}, objs[i]);

		bch_logger.log(LINE, " Worker ", thread_id," STOP *", LINE, "\n\n");
	}
	if (bch_logger.enable_logging_) { bch::Mutex.unlock(); }
}

// DONE
template<class bch_class>
void populateVectorOfMessagePolynomials(
		const unsigned char* str, 
		const threadZones& zones,
		std::vector<std::bitset<bch_class::k_>>& vector_of_message_polynomials)
{	
	std::bitset<bch_class::k_> temp_bitset;
	int it = zones.MPTG_end;
	int bit_pos = zones.bit_pos;
	for (int i = zones.MBTG_beginning; i < zones.MBTG_end; i++) {
		for (int j = 0; j < 8; ++j) {
			temp_bitset[bit_pos] = str[i] & (1 << j);
			bit_pos++;
			if (bit_pos == bch_class::k_) {
				vector_of_message_polynomials[it] ^= temp_bitset;
				it++;
				temp_bitset.reset();
				bit_pos = 0;
			}
		}
	}
	if (zones.bit_pos != 0) {
		vector_of_message_polynomials[it] ^= temp_bitset;
	}
}

template<class bch_class>
void populateUnsignedCharVectors(
		const std::vector<bchType>& objs,
		const threadZones& zones,
		const unsigned char* buffer)
{	
	int it = zones.MPTG_end;
	int bit_pos = zones.bit_pos;

	auto temp_poly = std::visit(
		[=](auto& obj) {
			using current_class = CVBC<decltype(obj.get())>;
			// becasuse in decodeBch we reversed it back to normal we need to move the data to the right 
			// so that the data bitset in union can read correct values
			obj->codeword_polynomials_.received.codeword >>= (current_class::n_ - current_class::k_);
			obj->codeword_polynomials_.decoded.codeword >>= (current_class::n_ - current_class::k_);
			return std::any(obj->codeword_polynomials_);
		}, objs[it]);
	auto polynomials = std::any_cast<polynomialData<bch_class::n_, bch_class::k_>>(temp_poly);

	for (int i = zones.MBTG_beginning; i < zones.MBTG_end; i++) {
		for (int j = 0; j < 8; ++j) {
			bch::modified_charstream[i] ^= polynomials.received.data[bit_pos] << j;
			bch::recovered_charstream[i] ^= polynomials.decoded.data[bit_pos] << j;
			bit_pos++;
			if (bit_pos == bch_class::k_) {
				it++;
				temp_poly = std::visit(
					[=](auto& obj) {
						using current_class = CVBC<decltype(obj.get())>;
						obj->codeword_polynomials_.received.codeword >>= (current_class::n_ - current_class::k_);
						obj->codeword_polynomials_.decoded.codeword >>= (current_class::n_ - current_class::k_);
						return std::any(obj->codeword_polynomials_);
					}, objs[it]);
				polynomials = std::any_cast<polynomialData<bch_class::n_, bch_class::k_>>(temp_poly);
				bit_pos = 0;
			}
		}
	}
}
// DONE
int parse_arguments(
		const int argc, 
		char* argv[])
{

	int fd = FAIL;
	int c;

	if (argc != 2 && argc != 7 && argc != 8) {
		std::cout << "Invalid number of arguments"<< std::endl;
		return FAIL;
	}

	while ((c = getopt (argc, argv, "i:c:p:vh")) != FAIL) {
		switch (c) {
			case 'h':
				printHelpMessage(argv[0]);
				exit(0);
			case 'i':
				bch::filename = optarg;
				fd = open(optarg, O_RDONLY);
				if (fd == FAIL) {
					std::cout << "Incorrect Image file name" << std::endl;
					return FAIL;
				}
				break;
			case 'c':
				if (atoi(optarg) < 0 || atoi(optarg) > 3) {
					return FAIL;
				}
				bch::code_type = static_cast<codeType>(atoi(optarg));
				break;
			case 'p':
				bch::error_probability = atoi(optarg);
				if (bch::error_probability < 10 || bch::error_probability > 10000000) {
					std::cout << "Error probability should be between (10 and 10000000)" << std::endl;
					return FAIL;
				}
				break;
			case 'v':
				bch_logger.enable_logging_ = true;
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

std::string bchFirstInit()
{
	switch (bch::code_type) {
		case (BCH6351):
			// initialize temp obj so that we can use visit for variant in later function
			bch::BCH_objects.push_back(std::move(std::make_unique<Bch6351>(0b0)));
			return "BCH6351.bmp";
		case (BCH6345):
			bch::BCH_objects.push_back(std::move(std::make_unique<Bch6345>(0b0)));
			return "BCH6345.bmp";
		case (BCH4836):
			bch::BCH_objects.push_back(std::move(std::make_unique<Bch4836>(0b0)));
			return "BCH4836.bmp";
		case (BCH4830):
			bch::BCH_objects.push_back(std::move(std::make_unique<Bch4830>(0b0)));
			return "BCH4830.bmp";
		default:
			return NULL;
	}
}

template<class bch_class>
void initializeBchMathStruct() 
{	
	bch::bch_math<bch_class::n_> = new mathHelper<bch_class::n_>;
	readPrimitivePolynomial<bch_class>(bch::bch_math<bch_class::n_>);
	generateGaloisField<bch_class>(bch::bch_math<bch_class::n_>); 
	generateGeneratorPolynomial<bch_class>(bch::bch_math<bch_class::n_>);
}

int main(
		int argc, 
		char* argv[])
{
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

	std::string file_suffix = bchFirstInit();

	std::string image_with_errors_path = bch::filename.substr(0, bch::filename.length()-4).append("_with_errors_" + file_suffix);
	std::string image_fixed_path = bch::filename.substr(0, bch::filename.length()-4).append("_fixed_" + file_suffix);

	remove(image_with_errors_path.c_str());
	remove(image_fixed_path.c_str());

	std::ofstream image_with_errors (image_with_errors_path.c_str(), std::ios::out | std::ios::app | std::ios::binary);
	std::ofstream image_fixed (image_fixed_path.c_str(), std::ios::out | std::ios::app | std::ios::binary);

	// resize vectors:
	const ssize_t NUMBER_OF_MESSAGE_POLYNOMIALS = std::visit([=](const auto& obj){
		using current_class = CVBC<decltype(obj.get())>;

		// do a check to see if there are any std::left over bites and if there are we need to add another message polynomial
		ssize_t NUMBER_OF_MESSAGE_POLYNOMIALS = ((FILE_BYTE_SIZE*8)%current_class::k_)? 
					(FILE_BYTE_SIZE*8)/current_class::k_ + 1 :
					(FILE_BYTE_SIZE*8)/current_class::k_;

		current_class::vector_of_message_polynomials.resize(NUMBER_OF_MESSAGE_POLYNOMIALS);

		return NUMBER_OF_MESSAGE_POLYNOMIALS;
	}, bch::BCH_objects[0]);
	bch::BCH_objects.resize(NUMBER_OF_MESSAGE_POLYNOMIALS);
	bch::recovered_charstream.resize(FILE_BYTE_SIZE);
	bch::modified_charstream.resize(FILE_BYTE_SIZE);

	std::vector<std::thread> threads;
	const int NUM_THREADS = std::thread::hardware_concurrency();
	std::cout << "number of detected threads: " << NUM_THREADS << std::endl;

	// MESSAGE_POLYNOMIALS_THREAD_GROUP --> count of group of message polynomials that will be used by a single thread
	const int MESSAGE_POLYNOMIALS_THREAD_GROUP = (NUMBER_OF_MESSAGE_POLYNOMIALS%NUM_THREADS)? 
													NUMBER_OF_MESSAGE_POLYNOMIALS/NUM_THREADS :
													NUMBER_OF_MESSAGE_POLYNOMIALS/NUM_THREADS - 1;

	// MESSAGE_BYTES_THREAD_GROUP --> count of group of message bytes that will be used by a single thread
	const int MESSAGE_BYTES_THREAD_GROUP = FILE_BYTE_SIZE / NUM_THREADS;

	std::cout << "Parsing image file..." << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	int old_bit_pos{0};
	int overlaping_message_polynomial{0};

	threadZones zones{};

	for (int thread_id = 0; thread_id < NUM_THREADS; thread_id++) {
		zones.MBTG_beginning = MESSAGE_BYTES_THREAD_GROUP*thread_id;
		zones.MBTG_end = MESSAGE_BYTES_THREAD_GROUP*(thread_id + 1);

		std::visit([&](const auto& obj){
			using current_class = CVBC<decltype(obj.get())>;
			zones.bit_pos = (zones.MBTG_beginning*8) % current_class::k_;
		}, bch::BCH_objects[0]);

		if (thread_id == NUM_THREADS-1) {
			zones.MBTG_end = FILE_BYTE_SIZE;
		}

		if (zones.bit_pos < old_bit_pos) {
			overlaping_message_polynomial++;
		}

		zones.MPTG_end = MESSAGE_POLYNOMIALS_THREAD_GROUP*thread_id + overlaping_message_polynomial;

		old_bit_pos = zones.bit_pos;

		std::visit([&](const auto& obj){
			using current_class = CVBC<decltype(obj.get())>;
			threads.push_back(std::thread(
										populateVectorOfMessagePolynomials<current_class>, 
										buffer, 
										zones,
										std::ref(current_class::vector_of_message_polynomials)));
		}, bch::BCH_objects[0]);
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

	std::visit([](const auto& obj){
		using current_class = CVBC<decltype(obj.get())>;
		initializeBchMathStruct<current_class>();
	}, bch::BCH_objects[0]);

	threads.clear();

	bch_logger.log("Please be patient, starting coding and decoding process...\n");

	for (int thread_id = 0; thread_id < NUM_THREADS; thread_id++) {
		zones.MPTG_beginning = MESSAGE_POLYNOMIALS_THREAD_GROUP*thread_id;
		zones.MPTG_end = MESSAGE_POLYNOMIALS_THREAD_GROUP*(thread_id + 1);

		if (thread_id == NUM_THREADS-1) {
			zones.MPTG_end = NUMBER_OF_MESSAGE_POLYNOMIALS;
		}
		
		threads.push_back(std::thread(beginMainProcess, std::ref(bch::BCH_objects), thread_id, zones));
		
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
		zones.MBTG_beginning = MESSAGE_BYTES_THREAD_GROUP*thread_id;
		zones.MBTG_end = MESSAGE_BYTES_THREAD_GROUP*(thread_id + 1);

		if (thread_id == NUM_THREADS-1) {
			zones.MBTG_end = FILE_BYTE_SIZE;
		}

		std::visit([&](const auto& obj){
			using current_class = CVBC<decltype(obj.get())>;
			zones.bit_pos = (zones.MBTG_beginning*8) % current_class::k_;
		}, bch::BCH_objects[0]);

		if (zones.bit_pos < old_bit_pos) { overlaping_message_polynomial++; }

		zones.MPTG_end = MESSAGE_POLYNOMIALS_THREAD_GROUP*thread_id + overlaping_message_polynomial;

		old_bit_pos = zones.bit_pos;

		std::visit([&](const auto& obj) -> void {
			using current_class = CVBC<decltype(obj.get())>;
			threads.push_back(std::thread(populateUnsignedCharVectors<current_class>,  std::ref(bch::BCH_objects), zones, buffer));
		}, bch::BCH_objects[0]);

	}

	for (auto& thread : threads) {
		thread.join();
	}

	stop = std::chrono::high_resolution_clock::now();
	dur = stop - start;
	duration = std::chrono::duration_cast<std::chrono::duration<float>>(dur);

	std::cout << "Converting modified and recovered data from bitsets to bytes done and it took: ";
	std::cout << std::to_string(duration.count()).substr(0, std::to_string(duration.count()).length()-3) << " seconds" << std::endl;

	// write bytes to new files and get diff
	for (int i = 0; i < HEADER_BYTES; i++) {
		image_with_errors << buffer[i];
		image_fixed << buffer[i];
	}

	for (ssize_t i = HEADER_BYTES; i < FILE_BYTE_SIZE; i++) {
		image_with_errors << bch::modified_charstream[i];
		image_fixed << bch::recovered_charstream[i];
	}

	size_t difference_count = 0;
	// TODO: maybe ignore first HEADER_BYTES so that it is a fair comparison
	for (int i = 0; i < NUMBER_OF_MESSAGE_POLYNOMIALS; i++) {
		std::visit([&](const auto& obj){
			using current_class = CVBC<decltype(obj.get())>;
			difference_count += (current_class::vector_of_message_polynomials[i] ^ 
									obj->codeword_polynomials_.decoded.data).count();
		}, bch::BCH_objects[i]);
	}

	auto table_row_separator = []()->void { std::cout << std::setw(42) << std::setfill('-') << "+" << std::setw(20) << "+" << std::setw(1) << "+" << std::setfill (' ') << std::endl; };
	auto table_items = [](const std::string &lhs, const std::string &rhs)->void { std::cout << "| " << std::setw(40) << lhs << std::setw(20) << "| " + rhs << "|" << std::endl; };

	std::visit([&](const auto& obj){
		using current_class = CVBC<decltype(obj.get())>;
		std::cout << std::endl << std::left;
		table_row_separator();
		table_items("Code used ", "(" + std::to_string(current_class::n_) + "," + std::to_string(current_class::k_) + "," + std::to_string(current_class::t_*2+1 ) + ")");
		table_row_separator();
		table_items("Original image used ", bch::filename);
		table_row_separator();
		table_items("Given error probability ", std::to_string(1/static_cast<float>(bch::error_probability)));
		table_row_separator();
		table_items("Real error probability ", std::to_string(1/static_cast<float>((NUMBER_OF_MESSAGE_POLYNOMIALS*current_class::n_)/counters.introduced_errors_count)));
		table_row_separator();
		table_items("Number of all data bits ", std::to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * current_class::k_));
		table_row_separator();
		table_items("Number of all data + redundant bits ", std::to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * current_class::n_));
		table_row_separator();
		table_items("Number of all introduced errors ", std::to_string(counters.introduced_errors_count));
		table_row_separator();
		table_items("Number of all message polynomials ", std::to_string(NUMBER_OF_MESSAGE_POLYNOMIALS));
		table_row_separator();
		table_items("Number of successful decodings", std::to_string(counters.success_count));
		table_row_separator();
		table_items("Number of decoding errors ", std::to_string(counters.failure_count));
		table_row_separator();
		table_items("Number of uncaught decoding errors ", std::to_string(counters.uncaught_errors_count));
		table_row_separator();
		table_items("Number of over t errors in codewords ", std::to_string(counters.big_errors_count));
		table_row_separator();
		table_items("Final data bits difference ", std::to_string(difference_count));
		table_row_separator();
		table_items("Encoding and decoding time ", std::to_string(main_duration.count()).substr(0, std::to_string(main_duration.count()).length()-3) + " seconds ");
		table_row_separator();
	}, bch::BCH_objects[0]);

	std::cout << std::endl << "Results have been written to BCH_logs.txt file" << std::endl;
	std::cout << std::endl << "To view made images on linux use \"feh -F -Z --force-aliasing -d " << bch::filename << " " << image_with_errors_path.c_str() << " " << image_fixed_path.c_str() << "\"" << std::endl << std::endl;

	std::ofstream BCH_logs;
	BCH_logs.open ("BCH_logs.txt", std::ios::out | std::ios::app);
	auto BCH_logs_items = [&](const std::string &lhs, const std::string &rhs)->void { BCH_logs << std::setw(37) << lhs << std::setw(20) << "| " + rhs << std::endl; };

	std::visit([&](const auto& obj){
		using current_class = CVBC<decltype(obj.get())>;
		BCH_logs << std::endl << std::left;
		BCH_logs_items("Code used ", "(" + std::to_string(current_class::n_) + "," + std::to_string(current_class::k_) + "," + std::to_string(current_class::t_*2+1 ) + ")");
		BCH_logs_items("Original image used ", bch::filename);
		BCH_logs_items("Given error probability ", std::to_string(1/static_cast<float>(bch::error_probability)));
		BCH_logs_items("Real error probability ", std::to_string(1/static_cast<float>((NUMBER_OF_MESSAGE_POLYNOMIALS * current_class::n_)/counters.introduced_errors_count)));
		BCH_logs_items("Number of all data bits ", std::to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * current_class::k_));
		BCH_logs_items("Number of all data + redundant bits ", std::to_string(NUMBER_OF_MESSAGE_POLYNOMIALS * current_class::n_));
		BCH_logs_items("Number of all generated errors ", std::to_string(counters.introduced_errors_count));
		BCH_logs_items("Number of all message polynomials ", std::to_string(NUMBER_OF_MESSAGE_POLYNOMIALS));
		BCH_logs_items("Number of successful decodings ", std::to_string(counters.success_count));
		BCH_logs_items("Number of decoding errors ", std::to_string(counters.failure_count));
		BCH_logs_items("Number of uncaught decoding errors", std::to_string(counters.uncaught_errors_count));
		BCH_logs_items("Number of over errors in codewords", std::to_string(counters.big_errors_count));
		BCH_logs_items("Final data bits difference ", std::to_string(difference_count));
		BCH_logs_items("Encoding and decoding time ", std::to_string(main_duration.count()).substr(0, std::to_string(main_duration.count()).length()-3) + " seconds ");
	}, bch::BCH_objects[0]);

	munmap(buffer, FILE_BYTE_SIZE);
	close(fd);

	image_with_errors.close();
	image_fixed.close();
	BCH_logs.close();

	return 0;
}
