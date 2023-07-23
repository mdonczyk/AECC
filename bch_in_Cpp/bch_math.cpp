#include "bch_math.hpp"

int multiplyIntPolynomials(
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