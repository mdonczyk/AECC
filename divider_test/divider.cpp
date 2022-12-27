#include <bitset>
#include <iostream>
#include <utility>

// Define the number of bits in the dividend and divisor.
constexpr int kNumBits = 63;

// Define a helper function to find the index of the most significant bit in a bitset.
int MSB(std::bitset<kNumBits> bits) {
  for (int i = kNumBits - 1; i >= 0; i--) {
    if (bits[i]) {
      return i;
    }
  }
  return -1;
}

// Define the function to perform the division.
std::pair<std::bitset<kNumBits>, std::bitset<kNumBits>> bitset_divider(
    std::bitset<kNumBits> dividend, std::bitset<kNumBits> divisor) {
  std::bitset<kNumBits> quotient;
  std::bitset<kNumBits> remainder = dividend;
  while (remainder.any() && (MSB(remainder) >= MSB(divisor))) {
    // Shift the divisor to align with the most significant bit of the remainder.
    int shift = MSB(remainder) - MSB(divisor);
    divisor <<= shift;
    // Update the quotient and remainder.
    quotient.set(shift);
    remainder ^= divisor;
    // Shift the divisor back to its original position.
    divisor >>= shift;
  }
  return {quotient, remainder};
}



int main() {
  // Define the dividend and divisor as bitsets.
  std::bitset<kNumBits> dividend(std::string("001000011111100000101001111100000000"));
  std::bitset<kNumBits> divisor(std::string("1001110010101"));

  // Perform the division and print the result.
  auto [quotient, remainder] = bitset_divider(dividend, divisor);
  std::cout << "Quotient: " << quotient << std::endl;
  std::cout << "Remainder: " << remainder << std::endl;

  return 0;
}
