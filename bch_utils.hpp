#ifndef BCH_UTILS_HPP
#define BCH_UTILS_HPP

#include <atomic>
#include <bit>
#include <bitset>
#include <iomanip>
#include <iostream>
#include <sys/stat.h>
#include <mutex>
#include <memory>
#include <fcntl.h>

#include "bch_logger.hpp"

#define GFB 63 // Galois Field Barrier = 2**m - 1 = 2**6 - 1

enum status {
    SUCCESS,
    FAIL    = -1
};

enum codeType {
    BCH6351,
    BCH6345,
    BCH4836,
    BCH4830
};

// global atomic counters:
struct globalCounters {
    std::atomic<size_t> success_count{0};
    std::atomic<size_t> failure_count{0};
    std::atomic<size_t> introduced_errors_count{0};
    std::atomic<size_t> big_errors_count{0};
    std::atomic<size_t> uncaught_errors_count{0};
};

struct threadZones {
	// message_bytes_thread_group beginning
	size_t MBTG_beginning;
	// message_bytes_thread_group end
	size_t MBTG_end;
	// message_polynomials_thread_group beginning
	size_t MPTG_beginning;
	// message_polynomials_thread_group end
	size_t MPTG_end;
	size_t bit_pos;
	size_t old_bit_pos;
	size_t overlaping_message_polynomial;
};

namespace bch {
    extern int error_probability;
    extern codeType code_type;
    extern std::string filename;
    extern std::mutex Mutex;
}

// TODO:  stylistic:
const std::string LINE(64/2 - 4, '*');
const std::string DASH_LINE(64/2, '-');

void printHelpMessage(
		const char *file_name);

int parse_arguments(
		const int argc,
		char* argv[]);

template <size_t N>
int MSB(
		const std::bitset <N> &Polynomial)
{
	return (N + (GFB-N) - static_cast<size_t>(std::countl_zero(Polynomial.to_ullong())));
}

template <size_t N>
void verbosePolynomial(
		const std::bitset <N> &Polynomial)
{
	int power = MSB(Polynomial);
	for (int i = power; i > 0; i--) {
		if (Polynomial[static_cast<size_t>(i)]) {
			if (i != power) {
				std::cout << " + " ;
			}
			if (i != 1) {
				std::cout << "x^" << i;
			} else {
				std::cout << "x";
			}
		}
	}
	std::cout << " + 1";
	std::cout << std::endl;
}

template<size_t N>
void reverseBitset(
		std::bitset <N> &Polynomial, 
		size_t Shift)
{
    for(size_t i = 0; i < N/2; ++i) {
    	bool temp_bit = Polynomial[i];
    	Polynomial[i] = Polynomial[N-i-1];
    	Polynomial[N-i-1] = temp_bit;
    }
	Polynomial >>= Shift;
}

template<size_t N>
std::bitset <N> tempReverseBitset(
		std::bitset <N> Polynomial, 
		size_t Shift)
{
    for(size_t i = 0; i < N/2; ++i) {
    	bool temp_bit = Polynomial[i];
    	Polynomial[i] = Polynomial[N-i-1];
    	Polynomial[N-i-1] = temp_bit;
    }
	return (Polynomial >>= Shift);
}

template <size_t N>
std::pair<std::bitset <N>, std::bitset <N>> divideBitsetPolynomials(
		const std::bitset <N> &dividend, 
		const std::bitset <N> &divisor)
{
	std::bitset <N> quotient, remainder = dividend;
 	while (MSB(remainder) >= MSB(divisor)) {
		size_t shift = static_cast<size_t>(MSB(remainder) - MSB(divisor));
		remainder ^= divisor << shift;
		quotient.flip(shift); 
	}
	return {remainder, quotient};
}

#endif /* BCH_UTILS_HPP */
