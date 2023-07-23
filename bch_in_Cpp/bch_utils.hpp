// #ifndef BCH_UTILS_HPP
// #define BCH_UTILS_HPP

#pragma once

#include <iostream>
#include <bit>
#include <bitset>
#include <atomic>
#include <unistd.h>
#include <sys/stat.h>
#include <mutex>

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
    std::atomic<int> success_count{0};
    std::atomic<int> failure_count{0};
    std::atomic<int> introduced_errors_count{0};
    std::atomic<int> big_errors_count{0};
    std::atomic<int> uncaught_errors_count{0};
};


struct threadZones {
	// MESSAGE_BYTES_THREAD_GROUP beginning
	ssize_t MBTG_beginning;
	// MESSAGE_BYTES_THREAD_GROUP end
	ssize_t MBTG_end;
	// MESSAGE_POLYNOMIALS_THREAD_GROUP beginning
	ssize_t MPTG_beginning;
	// MESSAGE_POLYNOMIALS_THREAD_GROUP end
	ssize_t MPTG_end;
	int bit_pos;
};


namespace bch {
    extern int error_probability;
    extern codeType code_type;
    extern std::string filename;
    extern std::mutex Mutex;
}

// TODO:  stylistic:
const std::string LINE(50/2 - 4, '*');
const std::string DASH_LINE(50/2, '-');

template <size_t N>
int MSB(
		const std::bitset <N> &Polynomial)
{
	return (N + (GFB-N) - std::countl_zero(Polynomial.to_ullong()));
}

template <size_t N>
void verbosePolynomial(
		const std::bitset <N> &Polynomial)
{
	int power = MSB(Polynomial);
	for (int i=power; i>=0; i--) {
		if (Polynomial[i]) {
			if (i != power) {
				std::cout << " + " ;
			}
			if (i!=0) {
				if (i == 1) {
					std::cout << "N";
				} else {
				std::cout << "N^" << i;
				}
			} else {
				std::cout << "1";
			}
		}
	}
	std::cout << std::endl;
}

template<size_t N>
void reverseBitset(
		std::bitset <N> &Polynomial, 
		int Shift)
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
		int Shift)
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
		int shift = MSB(remainder) - MSB(divisor);
		remainder ^= divisor << shift;
		quotient.flip(shift); 
	}
	return {remainder, quotient};
}

void printHelpMessage(const char *file_name);

// #endif /* BCH_UTILS_HPP */