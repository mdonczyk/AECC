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
#include <ranges>
#include <thread>
#include <time.h>
#include <unordered_map>
#include <variant>
#include <vector>
#include <sys/mman.h>
#include <sys/stat.h>

#include "bch_utils.hpp"
#include "bch_math.hpp"
#include "bch_logger.hpp"

#define RESERVED_BYTES 71

template <size_t N, size_t K>
struct polynomialData{
    union polynomial {
        std::bitset<N> codeword;
        std::bitset<K> data;
    };
    polynomial encoded{};
    polynomial received{};
    polynomial decoded{};
};

std::string bchFirstInit();

size_t resizeMainVectors(
        const size_t file_byte_size);

void divideImageBytesToBitsets(
		const char* buffer,
		const size_t message_bytes_thread_group,
		const size_t message_polynomials_thread_group,
		const size_t num_threads,
		const size_t file_byte_size,
		std::vector<std::thread>& threads);

void mathStructInit();

void startMainProcess(
		const size_t number_of_message_polynomials,
		const size_t message_polynomials_thread_group,
		const size_t num_threads,
		std::vector<std::thread>& threads);

void divideImageBitsetsToBytes(
		const size_t message_bytes_thread_group,
		const size_t message_polynomials_thread_group,
		const size_t num_threads,
		const size_t file_byte_size,
		std::vector<std::thread>& threads);

size_t getMainDifferenceCount(
		const size_t number_of_message_polynomials);

globalCounters& getFullCounters();

void finalLogsAndCleanup(
		const size_t number_of_message_polynomials,
		std::string& image_with_errors_path,
		std::string& image_fixed_path,
		std::chrono::duration<float> main_duration);

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
    inline std::vector <bchType> BCH_objects;
    inline std::vector <char> decoded_charstream;
    inline std::vector <char> received_charstream;
};

#endif /* BCH_SIMULATOR_HPP */
