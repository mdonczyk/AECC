#ifndef BCH_SIMULATOR_HPP
#define BCH_SIMULATOR_HPP

#include <any>
#include <bit>
#include <bitset>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <functional>
#include <memory>
#include <random>
#include <thread>
#include <time.h>
#include <unordered_map>
#include <variant>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#define RESERVED_BYTES 54

template <size_t N, size_t K>
struct polynomialData{
    union enc {
        std::bitset<N> codeword;
        std::bitset<K> data;
    };
    union rec {
        std::bitset<N> codeword;
        std::bitset<K> data;
    };
    union dec {
        std::bitset<N> codeword;
        std::bitset<K> data;
    };
    enc  encoded{};
    rec  received{};
    dec  decoded{};
};

class Bch6351 {
    public:
        static constexpr ssize_t n_ = 63;
        static constexpr ssize_t k_ = 51;
        static constexpr ssize_t t_ = (n_ - k_) / 6;
        explicit Bch6351 (const std::bitset <k_>& data) {
            codeword_polynomials_.encoded.data = data;
        }
        static std::vector <std::bitset <k_>> vector_of_message_polynomials;
        polynomialData<n_, k_> codeword_polynomials_{};
};

class Bch6345 {
    public:
        static constexpr ssize_t n_ = 63;
        static constexpr ssize_t k_ = 45;
        static constexpr ssize_t t_ = (n_ - k_) / 6;
        explicit Bch6345(const std::bitset <k_>& data) {
            codeword_polynomials_.encoded.data = data;
        }
        static std::vector <std::bitset <k_>> vector_of_message_polynomials;
        polynomialData<n_, k_> codeword_polynomials_{};
};

class Bch4836 {
    public:
        static constexpr ssize_t n_ = 48;
        static constexpr ssize_t k_ = 36;
        static constexpr ssize_t t_ = (n_ - k_) / 6;
        explicit Bch4836(const std::bitset <k_>& data) {
            codeword_polynomials_.encoded.data = data;
        }
        static std::vector <std::bitset <k_>> vector_of_message_polynomials;
        polynomialData<n_, k_> codeword_polynomials_{};
};

class Bch4830{
    public:
        static constexpr ssize_t n_ = 48;
        static constexpr ssize_t k_ = 30;
        static constexpr ssize_t t_ = (n_ - k_) / 6;
        explicit Bch4830(const std::bitset <k_>& data) {
            codeword_polynomials_.encoded.data = data;
        }
        static std::vector <std::bitset <k_>> vector_of_message_polynomials;
        polynomialData<n_, k_> codeword_polynomials_{};
};

using bchType = std::variant<std::unique_ptr<Bch6351>, std::unique_ptr<Bch6345>, 
                             std::unique_ptr<Bch4836>, std::unique_ptr<Bch4830>>;

namespace bch {
    // the big dawg
    std::vector <bchType> BCH_objects;
    std::vector <unsigned char> decoded_charstream;
    std::vector <unsigned char> received_charstream;
};

#endif /* BCH_SIMULATOR_HPP */
