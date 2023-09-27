#ifndef BCH_MATH_HPP
#define BCH_MATH_HPP

#include <bitset>
#include <memory>
#include <vector>
#include <iostream>
#include <set>

#include "bch_utils.hpp"

int multiplyIntPolynomials(
		int mulitplicand,
		int multiplicator);

template <size_t N>
struct mathHelper {
    int m;
    std::bitset <N> primitive_polynomial_bitset {0b1011011};
    int primitive_polynomial_int;
    int alpha_to[GFB], index_of[GFB];
    std::bitset <N> generator_polynomial_bitset;
    std::vector <int> zeros, g, errpos;
    std::vector <std::vector <int>> zeros_cosets;
};

namespace bch {
    template <size_t N>
    std::unique_ptr<mathHelper<N>> math_struct;
}

template <class bch_class>
void readPrimitivePolynomial(
		std::unique_ptr<mathHelper<bch_class::n_>>& math_struct)
{
	math_struct->primitive_polynomial_int = static_cast<int>(math_struct->primitive_polynomial_bitset.to_ulong());
	std::cout << "Primitive polynomial:" << std::endl << "p(x) = ";
	verbosePolynomial(math_struct->primitive_polynomial_bitset);
}

template<class bch_class>
void generateGaloisField(
		std::unique_ptr<mathHelper<bch_class::n_>>& math_struct)
{
	math_struct->index_of[0] = -1;
	math_struct->m = static_cast<int>(MSB(math_struct->primitive_polynomial_bitset));
	for (int i = 0; i < GFB; i++) {
		if (i < math_struct->m) {
			math_struct->alpha_to[i] = 1 << i;
			math_struct->index_of[math_struct->alpha_to[i]] = i;
		} else {
			if (math_struct->alpha_to[i - 1] >= 32) {
				math_struct->alpha_to[i] = (math_struct->alpha_to[i - 1] << 1) ^ 
										math_struct->primitive_polynomial_int;
			} else {
				math_struct->alpha_to[i] = math_struct->alpha_to[i - 1] << 1;
			}
			math_struct->index_of[math_struct->alpha_to[i]] = i;
		}
	}
}

template <class bch_class>
void generateGeneratorPolynomial(
		std::unique_ptr<mathHelper<bch_class::n_>>& math_struct)
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
		using index_t = std::vector<int>::size_type;
		for (index_t i = 1; i < static_cast<size_t>(math_struct->m); i++) {
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
			for (int root = 1; root <= 2*bch_class::t_; root++) {
				if (element == root) {
					root_found = true;
					break;
				}
			}
			if(root_found) {
				math_struct->zeros_cosets.push_back(coset);
				break; 
			}
		}
		if (math_struct->zeros_cosets.size() == bch_class::t_) {
			break;
		}
	}
	std::vector <int> min_polynomials;
	for (const auto& zero_coset : math_struct->zeros_cosets) {
		int product = math_struct->alpha_to[zero_coset[0]] ^ 2; // (ax + x)
		for (uint i = 1; i < zero_coset.size(); i++) {
			product = multiplyIntPolynomials(product, math_struct->alpha_to[zero_coset[i]] ^ 2);
		}
		product %= GFB+1;
		product ^= math_struct->primitive_polynomial_int;
		min_polynomials.push_back(product);
	}
	if (min_polynomials.empty()){
		exit(1);
	}
	int temp_poly = min_polynomials[0];
	for (size_t i=1; i<bch_class::t_; i++) {
		temp_poly = multiplyIntPolynomials(temp_poly, min_polynomials[i]);
	}
	math_struct->generator_polynomial_bitset = static_cast<unsigned>(temp_poly);
	reverseBitset(math_struct->generator_polynomial_bitset, bch_class::k_-1);
	std::cout << "This is a (" << bch_class::n_ << "," << bch_class::k_ << "," << bch_class::t_*2+1 << ") binary bch code" << std::endl;
	std::cout << "g(x) is " << math_struct->generator_polynomial_bitset.to_string().substr(bch_class::n_ - (bch_class::n_-bch_class::k_+1)) << std::endl;
}

#endif /* math_struct_HPP */
