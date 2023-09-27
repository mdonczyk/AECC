#include "bch_simulator.hpp"

using bchType = std::variant<std::unique_ptr<Bch6351>, std::unique_ptr<Bch6345>, 
                             std::unique_ptr<Bch4836>, std::unique_ptr<Bch4830>>;

template<typename T>
// (CURRENT VARIANT BCH CLASS) Alias for getting the actual type of the Bch code type class used in variant
using CVBC = std::remove_pointer_t<std::decay_t<T>>;

std::vector<std::bitset <Bch6351::k_>> Bch6351::vector_of_message_polynomials;
std::vector<std::bitset <Bch6345::k_>> Bch6345::vector_of_message_polynomials;
std::vector<std::bitset <Bch4836::k_>> Bch4836::vector_of_message_polynomials;
std::vector<std::bitset <Bch4830::k_>> Bch4830::vector_of_message_polynomials;

static globalCounters counters{};

std::string bchFirstInit()
{
	switch (bch::code_type) {
		case (BCH6351):
			// initialize temp obj so that we can use visit for variant in later function
			bch::BCH_objects.push_back(std::make_unique<Bch6351>(0b0));
			return "BCH6351.bmp";
		case (BCH6345):
			bch::BCH_objects.push_back(std::make_unique<Bch6345>(0b0));
			return "BCH6345.bmp";
		case (BCH4836):
			bch::BCH_objects.push_back(std::make_unique<Bch4836>(0b0));
			return "BCH4836.bmp";
		case (BCH4830):
			bch::BCH_objects.push_back(std::make_unique<Bch4830>(0b0));
			return "BCH4830.bmp";
		default:
			return NULL;
	}
}

size_t resizeMainVectors(const size_t file_byte_size) {
	const size_t number_of_message_polynomials = std::visit([=](const auto& obj){
		using current_class = CVBC<decltype(obj.get())>;

		// do a check to see if there are any left over bites and if
		// there are we need to add another message polynomial
		size_t NOMP = ((file_byte_size*8)%current_class::k_)? 
					(file_byte_size*8)/current_class::k_ + 1 :
					(file_byte_size*8)/current_class::k_;

		current_class::vector_of_message_polynomials.resize(NOMP);

		return NOMP;
	}, bch::BCH_objects[0]);

	bch::BCH_objects.resize(number_of_message_polynomials);
	bch::decoded_charstream.resize(file_byte_size);
	bch::received_charstream.resize(file_byte_size);

	return number_of_message_polynomials;
}

template<class bch_class>
void populateVectorOfMessagePolynomials(
		const char* buffer, 
		const threadZones& zones,
		std::vector<std::bitset<bch_class::k_>>& vector_of_message_polynomials)
{	
	size_t data_bitset_iterator = zones.MPTG_end;
	size_t bit_pos = zones.bit_pos;
	for (size_t i = zones.MBTG_beginning; i < zones.MBTG_end; i++) {
		for (int j = 0; j < 8; ++j) {
			vector_of_message_polynomials[data_bitset_iterator][bit_pos] = buffer[i] & (1 << j);
			bit_pos++;
			if (bit_pos == bch_class::k_) {
				data_bitset_iterator++;
				bit_pos = 0;
			}
		}
	}
}

void divideImageBytesToBitsets(
		const char* buffer,
		const size_t message_bytes_thread_group,
		const size_t message_polynomials_thread_group,
		const size_t num_threads,
		const size_t file_byte_size,
		std::vector<std::thread>& threads)
{
	threadZones zones{};

	for (size_t thread_id = 0; thread_id < num_threads; thread_id++) {
		zones.MBTG_beginning = message_bytes_thread_group*thread_id;
		zones.MBTG_end = message_bytes_thread_group*(thread_id + 1);

		std::visit([&](const auto& obj){
			using current_class = CVBC<decltype(obj.get())>;
			zones.bit_pos = (zones.MBTG_beginning*8) % current_class::k_;
		}, bch::BCH_objects[0]);

		if (thread_id == num_threads-1) {
			zones.MBTG_end = file_byte_size;
		}

		if (zones.bit_pos < zones.old_bit_pos) {
			zones.overlaping_message_polynomial++;
		}

		zones.MPTG_end = message_polynomials_thread_group*thread_id + zones.overlaping_message_polynomial;

		zones.old_bit_pos = zones.bit_pos;

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
	
	threads.clear();
}

template<class bch_class>
void encodeBch(
		polynomialData<bch_class::n_, bch_class::k_>& polynomials, 
		const std::unique_ptr<mathHelper<bch_class::n_>>& math_struct)
{	
	polynomials.encoded.codeword <<= (bch_class::n_-bch_class::k_);
	std::bitset<bch_class::n_> rb = divideBitsetPolynomials(polynomials.encoded.codeword, math_struct->generator_polynomial_bitset).first;
	polynomials.encoded.codeword ^= rb;
	reverseBitset(polynomials.encoded.codeword, 0);
	bch_logger::log("\n");
}

template<class bch_class>
std::vector <int> calculateSyndromes(
		const polynomialData<bch_class::n_, bch_class::k_>& polynomials, 
		const std::unique_ptr<mathHelper<bch_class::n_>>& math_struct, 
		bool &errors_in_codeword)
{
	// for (auto& syndrome : syndromes | std::ranges::views::drop(1)) {
	std::vector <int> syndromes(2*bch_class::t_+1);
	for (size_t i = 1; i <= 2*bch_class::t_; i++) {
		syndromes[i] = 0;
		for (size_t j = 0; j < bch_class::n_; j++) {
			if (polynomials.received.codeword[j] != 0) {
				syndromes[i] ^= math_struct->alpha_to[(i*j) % GFB];
			}
		}
		syndromes[i] = math_struct->index_of[syndromes[i]];
		// when there are no errors all syndromes should be -1
		if (syndromes[i] != -1) {
			errors_in_codeword = true;
		}
	}

	bch_logger::log("syndromes= ( ");
		for (std::vector<int>::size_type i = 1; i <= 2*bch_class::t_; i++) {
			bch_logger::log(syndromes[i], " ");
		}
	bch_logger::log(")\n");

	return syndromes;
}

template<class bch_class>
void introduceErrors(
		polynomialData<bch_class::n_, bch_class::k_>& polynomials) 
{
	polynomials.received.codeword = polynomials.encoded.codeword;
	std::random_device rd;
	std::mt19937 gen (rd());
	std::uniform_int_distribution<> d (0, bch::error_probability);
	for (size_t i = 0; i < bch_class::n_; i++) {
		int randnum = d(gen);
		if (randnum == 0) {
			counters.introduced_errors_count ++;
			polynomials.received.codeword.flip(i);
		}
	}
}

template<class bch_class>
void compareAndPrintCodewords(
		const polynomialData<bch_class::n_, bch_class::k_>& polynomials)
{	
	bch_logger::log("c(x) = ", tempReverseBitset(polynomials.encoded.codeword, 0), "\n");
	bch_logger::log("r(x) = ", tempReverseBitset(polynomials.received.codeword, 0), "\n");

	std::bitset <bch_class::n_> bit_difference = polynomials.encoded.codeword ^ polynomials.received.codeword;
	reverseBitset(bit_difference, 0);

	if (bit_difference.count() != 0) {
		if (bch_logger::enable_logging) {
			bch_logger::log("       ", (bit_difference.to_string(' ','^')), "\n");
			bch_logger::log("Position of ", bit_difference.count());
			bit_difference.count() == 1? bch_logger::log(" error "): bch_logger::log(" errors ");
			bch_logger::log("in the received codeword at ^.\n");
		}
		if (bit_difference.count() > static_cast<unsigned>(bch_class::t_)) {
			counters.big_errors_count++;
		}
	} else {
		bch_logger::log("No errors in the received codeword.\n");
	}
}

template<class bch_class>
status decodeBch(
		polynomialData<bch_class::n_, bch_class::k_>& polynomials,
		const std::unique_ptr<mathHelper<bch_class::n_>>& math_struct) 
{
	/*
	 Simon Rockliff's implementation of Berlekamp's algorithm.
	 Assume we have received bits in polynomials.received.codeword[i], i=0..(bch_class::n_-1).
	
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
	std::vector <int> syndromes = calculateSyndromes<bch_class>(polynomials, math_struct, errors_in_codeword);
	polynomials.decoded.codeword = polynomials.received.codeword;
	
	// reverse bitsets back to normal and shift data to the right:
	reverseBitset(polynomials.encoded.codeword, 0);
	reverseBitset(polynomials.received.codeword, 0);
	polynomials.encoded.codeword >>= bch_class::n_ - bch_class::k_;
	polynomials.received.codeword >>= bch_class::n_ - bch_class::k_;

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
					elp[u][i] = math_struct->index_of[elp[u][i]];
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
						elp[u + 1][i + u - q] = math_struct->alpha_to[(d[u] + 
							GFB - d[q] + elp[q][i]) % GFB];
					}
				}
				for (int i = 0; i <= l[u]; i++) {
					elp[u + 1][i] ^= elp[u][i];
					elp[u][i] = math_struct->index_of[elp[u][i]];
				}
			}
			u_lu[u + 1] = u - l[u + 1];
			// form (u+1)th discrepancy
			if (u < 2*bch_class::t_) {
			// no discrepancy computed on last iteration
				using index_t = std::vector<int>::size_type;
				index_t u_t = static_cast<index_t>(u);
				if (syndromes[u_t + 1] != -1) {
					d[u + 1] = math_struct->alpha_to[syndromes[u_t + 1]];
				} else {
				d[u + 1] = 0;
				}
				for (int i = 1; i <= l[u + 1]; i++) {
					if ((syndromes[u_t + 1 - static_cast<index_t>(i)] != -1) && (elp[u + 1][i] != 0)) {
						d[u + 1] ^= math_struct->alpha_to[(syndromes[u_t + 1 - static_cast<index_t>(i)] + 
							math_struct->index_of[elp[u + 1][i]]) % GFB];
					}
				}
				// put d[u+1] into index form
				d[u + 1] = math_struct->index_of[d[u + 1]];
			}
		} while ((u < 2*bch_class::t_) && (l[u + 1] <= bch_class::t_));
		u++;
		if (l[u] <= bch_class::t_) {// Can correct errors
			// put elp into index form
			bch_logger::log("Sigma(x) = ");
			for (int i = 0; i <= l[u]; i++) {
				elp[u][i] = math_struct->index_of[elp[u][i]];
				bch_logger::log(elp[u][i], " ");
			}
			// Find roots of the error location polynomial:
			// Chien search
			bch_logger::log("\nRoots: ");
			for (int i = 1; i <= GFB; i++) {
				int q = 1;
				for (int j = 1; j <= l[u]; j++) {
					elp[u][j] = (elp[u][j] + j) % GFB;
					q ^= math_struct->alpha_to[elp[u][j]];
				}
				// Store error location number indices
				if (!q) {
					error_locations.push_back(i  - 1 -(GFB-bch_class::n_));
					bch_logger::log(error_locations.back(), " ");
				}
			}
			bch_logger::log("\n");
			if (error_locations.size() == static_cast<unsigned>(l[u])) {
				for (auto const &error_location : error_locations) {
					size_t err_loc = static_cast<size_t>((bch_class::n_-1 - error_location + bch_class::n_) % bch_class::n_);
					if (err_loc < bch_class::n_) {
						polynomials.decoded.codeword.flip(err_loc);
					} else {
						bch_logger::log("\nIncomplete decoding: errors detected (err_loc out of bound: err_loc = ", err_loc, ")\n");
						goto decode_error;
					}
				}
			} else {
				bch_logger::log("\nIncomplete decoding: errors detected (elp has degree > bch_class::t_)\n");
				goto decode_error;
			}
		} else {
			bch_logger::log("\nIncomplete decoding: errors detected (l[u] > bch_class::t_: l[u] = ", l[u], ")\n");
			goto decode_error;
		}
	} else {
		bch_logger::log("No errors found\n");
	}
	reverseBitset(polynomials.decoded.codeword, 0);
	polynomials.decoded.codeword >>= bch_class::n_ - bch_class::k_;
	return SUCCESS;

	decode_error:
		reverseBitset(polynomials.decoded.codeword, 0);
		polynomials.decoded.codeword >>= bch_class::n_ - bch_class::k_;
		return FAIL;
}

template<class bch_class>
void compareAndPrintData(
		const polynomialData<bch_class::n_, bch_class::k_>& polynomials)
{
	bch_logger::log(DASH_LINE, "Results", DASH_LINE,"\n");

	bch_logger::log("Original Data  = ", polynomials.encoded.data, "\n");
	bch_logger::log("Recovered Data = ", polynomials.decoded.data, "\n");

	std::bitset <bch_class::k_> bit_difference = polynomials.encoded.data ^ polynomials.decoded.data;

	if (bit_difference.count() != 0) {
		bch_logger::log("                 ", (bit_difference.to_string(' ','^')), "\n");
		bch_logger::log("Position of ", bit_difference.count()," message decoding");
		bit_difference.count() == 1? bch_logger::log(" error "): bch_logger::log(" errors ");
		bch_logger::log("at ^.\n");
		counters.uncaught_errors_count++;
	} else {
		bch_logger::log("Successful decoding.\n");
	}
}

void mainProcess(
		std::vector<bchType>& objs,
		const int thread_id, 
		const threadZones& zones)
{	
	if (bch_logger::enable_logging) { bch::Mutex.lock(); }
	for (size_t i = zones.MPTG_beginning; i < zones.MPTG_end; i++) {
		bch_logger::log(LINE, " Worker ", thread_id," START ", LINE);

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

			encodeBch<current_class>(obj->codeword_polynomials_, bch::math_struct<current_class::n_>);

			introduceErrors<current_class>(obj->codeword_polynomials_);

			compareAndPrintCodewords<current_class>(obj->codeword_polynomials_);

			if (status st = decodeBch<current_class>(obj->codeword_polynomials_, bch::math_struct<current_class::n_>); st == SUCCESS) {
				compareAndPrintData<current_class>(obj->codeword_polynomials_);
				counters.success_count++;
			} else {
				counters.failure_count++;
			}
		}, objs[i]);

		bch_logger::log(LINE, " Worker ", thread_id," STOP *", LINE, "\n\n");
	}
	if (bch_logger::enable_logging) { bch::Mutex.unlock(); }
}


void startMainProcess(
		const size_t number_of_message_polynomials,
		const size_t message_polynomials_thread_group,
		const size_t num_threads,
		std::vector<std::thread>& threads)
{
	threadZones zones = {};

	for (size_t thread_id = 0; thread_id < num_threads; thread_id++) {
		zones.MPTG_beginning = message_polynomials_thread_group*thread_id;
		zones.MPTG_end = message_polynomials_thread_group*(thread_id + 1);

		if (thread_id == num_threads-1) {
			zones.MPTG_end = number_of_message_polynomials;
		}
		
		threads.push_back(std::thread(
									mainProcess,
									std::ref(bch::BCH_objects),
									thread_id,
									zones));
	}

	for (auto& thread : threads) {
		thread.join();
	}

	threads.clear();
}

template<class bch_class>
void populateCharVectors(
		const std::vector<bchType>& objs,
		const threadZones& zones)
{	
	size_t data_bitset_iterator = zones.MPTG_end;
	size_t bit_pos = zones.bit_pos;

	auto temp_poly = std::visit(
		[=](auto& obj) {
			return std::any(obj->codeword_polynomials_);
		}, objs[data_bitset_iterator]);
	auto polynomials = std::any_cast<polynomialData<bch_class::n_, bch_class::k_>>(temp_poly);

	for (size_t i = zones.MBTG_beginning; i < zones.MBTG_end; i++) {
		for (int j = 0; j < 8; j++) {
			bch::received_charstream[i] ^= static_cast<char>(polynomials.received.data[bit_pos] << j);
			bch::decoded_charstream[i] ^= static_cast<char>(polynomials.decoded.data[bit_pos] << j);
			bit_pos++;
			if (bit_pos == bch_class::k_) {
				data_bitset_iterator++;
				temp_poly = std::visit(
					[=](auto& obj) {
						return std::any(obj->codeword_polynomials_);
					}, objs[data_bitset_iterator]);
				polynomials = std::any_cast<polynomialData<bch_class::n_, bch_class::k_>>(temp_poly);
				bit_pos = 0;
			}
		}
	}
}

void divideImageBitsetsToBytes(
		const size_t message_bytes_thread_group,
		const size_t message_polynomials_thread_group,
		const size_t num_threads,
		const size_t file_byte_size,
		std::vector<std::thread>& threads)
{
    threadZones zones = {};

	for (size_t thread_id = 0; thread_id < num_threads; thread_id++) {
		zones.MBTG_beginning = message_bytes_thread_group*thread_id;
		zones.MBTG_end = message_bytes_thread_group*(thread_id + 1);

        std::visit([&](const auto& obj){
			using current_class = CVBC<decltype(obj.get())>;
			zones.bit_pos = (zones.MBTG_beginning*8) % current_class::k_;
		}, bch::BCH_objects[0]);

		if (thread_id == num_threads-1) {
			zones.MBTG_end = file_byte_size;
		}

		if (zones.bit_pos < zones.old_bit_pos) { 
            zones.overlaping_message_polynomial++;
        }

		zones.MPTG_end = message_polynomials_thread_group*thread_id + zones.overlaping_message_polynomial;

		zones.old_bit_pos = zones.bit_pos;

		std::visit([&](const auto& obj) -> void {
			using current_class = CVBC<decltype(obj.get())>;
			threads.push_back(std::thread(
                                        populateCharVectors<current_class>,
                                        std::ref(bch::BCH_objects),
                                        zones));
		}, bch::BCH_objects[0]);

	}

	for (auto& thread : threads) {
		thread.join();
	}

    threads.clear();
}

template<class bch_class>
void initializeBchMathStruct() 
{	
	bch::math_struct<bch_class::n_> = std::make_unique<mathHelper<bch_class::n_>>();
	readPrimitivePolynomial<bch_class>(bch::math_struct<bch_class::n_>);
	generateGaloisField<bch_class>(bch::math_struct<bch_class::n_>); 
	generateGeneratorPolynomial<bch_class>(bch::math_struct<bch_class::n_>);
}

void mathStructInit() 
{
	std::visit([](const auto& obj){
		using current_class = CVBC<decltype(obj.get())>;
		initializeBchMathStruct<current_class>();
	}, bch::BCH_objects[0]);
}

void finalLogsAndCleanup(
		const size_t number_of_message_polynomials,
		std::string& image_with_errors_path,
		std::string& image_fixed_path,
		std::chrono::duration<float> main_duration)
{
	size_t difference_count{0};

	for (size_t i = 0; i < number_of_message_polynomials; i++) {
		std::visit([&](const auto& obj){
			using current_class = CVBC<decltype(obj.get())>;
			difference_count += (current_class::vector_of_message_polynomials[i] ^ 
									obj->codeword_polynomials_.decoded.data).count();
		}, bch::BCH_objects[i]);
	}

	const auto table_row_separator = []()->void {
		std::cout << std::setw(42) << std::setfill('-') << "+" << std::setw(20) << "+" << std::setw(1) << "+" << std::setfill (' ') << std::endl;
	};

	const auto table_items = [](const std::string &lhs, const std::string &rhs)->void {
		std::cout << "| " << std::setw(40) << lhs << std::setw(20) << "| " + rhs << "|" << std::endl;
	};

	
	std::visit([&](const auto& obj){
		using current_class = CVBC<decltype(obj.get())>;
		std::cout << std::endl << std::left;
		table_row_separator();
		table_items("Code used ", "( " + std::to_string(current_class::n_) + ", " + std::to_string(current_class::k_) + ", " + std::to_string(current_class::t_*2+1 ) + " )");
		table_row_separator();
		table_items("Original image used ", bch::filename);
		table_row_separator();
		table_items("Given error probability ", std::to_string(1/static_cast<float>(bch::error_probability)));
		table_row_separator();
		if (counters.introduced_errors_count != 0) {
			table_items("Real error probability ", std::to_string(1/static_cast<float>((number_of_message_polynomials*current_class::n_)/counters.introduced_errors_count)));
			table_row_separator();
		}
		table_items("Number of all data bits ", std::to_string(number_of_message_polynomials * current_class::k_));
		table_row_separator();
		table_items("Number of all data + redundant bits ", std::to_string(number_of_message_polynomials * current_class::n_));
		table_row_separator();
		table_items("Number of all introduced errors ", std::to_string(counters.introduced_errors_count));
		table_row_separator();
		table_items("Number of all message polynomials ", std::to_string(number_of_message_polynomials));
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

		std::cout << std::endl << "Results have been written to BCH_logs.txt file" << std::endl;
		std::cout << std::endl << "To view made images on linux use \"feh -F -Z --force-aliasing -d " << bch::filename + " " + image_with_errors_path + " " + image_fixed_path + "\"" << std::endl << std::endl;

		std::ofstream BCH_logs;
		BCH_logs.open ("BCH_logs.txt", std::ios::out | std::ios::app);
		auto BCH_logs_items = [&](const std::string &lhs, const std::string &rhs)->void { BCH_logs << std::setw(37) << lhs << std::setw(20) << "| " + rhs << std::endl; };

		BCH_logs << std::endl << std::left;
		BCH_logs_items("Code used ", "( " + std::to_string(current_class::n_) + ", " + std::to_string(current_class::k_) + ", " + std::to_string(current_class::t_*2+1 ) + " )");
		BCH_logs_items("Original image used ", bch::filename);
		BCH_logs_items("Given error probability ", std::to_string(1/static_cast<float>(bch::error_probability)));
		if (counters.introduced_errors_count != 0) {
			BCH_logs_items("Real error probability ", std::to_string(1/static_cast<float>((number_of_message_polynomials * current_class::n_)/counters.introduced_errors_count)));
		}
		BCH_logs_items("Number of all data bits ", std::to_string(number_of_message_polynomials * current_class::k_));
		BCH_logs_items("Number of all data + redundant bits ", std::to_string(number_of_message_polynomials * current_class::n_));
		BCH_logs_items("Number of all generated errors ", std::to_string(counters.introduced_errors_count));
		BCH_logs_items("Number of all message polynomials ", std::to_string(number_of_message_polynomials));
		BCH_logs_items("Number of successful decodings ", std::to_string(counters.success_count));
		BCH_logs_items("Number of decoding errors ", std::to_string(counters.failure_count));
		BCH_logs_items("Number of uncaught decoding errors", std::to_string(counters.uncaught_errors_count));
		BCH_logs_items("Number of over errors in codewords", std::to_string(counters.big_errors_count));
		BCH_logs_items("Final data bits difference ", std::to_string(difference_count));
		BCH_logs_items("Encoding and decoding time ", std::to_string(main_duration.count()).substr(0, std::to_string(main_duration.count()).length()-3) + " seconds ");

		BCH_logs.close();
	}, bch::BCH_objects[0]);
}
