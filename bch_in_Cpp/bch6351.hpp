#include <any>
#include <atomic>
#include <bit>
#include <bitset>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <functional>
#include <memory>
#include <mutex>
#include <random>
#include <set>
#include <thread>
#include <time.h>
#include <unordered_map>
#include <variant>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#define GFB 63 // Galois Field Barrier = 2**m - 1 = 2**6 - 1
#define HEADER_BYTES 30

// global atomic counters:
struct globalCounters {
    std::atomic<int> success_count{0};
    std::atomic<int> failure_count{0};
    std::atomic<int> introduced_errors_count{0};
    std::atomic<int> big_errors_count{0};
    std::atomic<int> uncaught_errors_count{0};
};


// TODO:  stylistic:
const std::string LINE(50/2 - 4, '*');
const std::string DASH_LINE(50/2, '-');

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

// forward declaration for mathHelper:
template <size_t N>
struct mathHelper;
template <size_t N>
void readPrimitivePolynomial(mathHelper<N>& bch_math);
template <size_t N>
void generateGaloisField(mathHelper<N>& bch_math);
template <size_t N>
void generateGeneratorPolynomial(mathHelper<N>& bch_math);

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

class Bch6351 {
    public:
        explicit Bch6351 (const std::bitset <51> &data) {
            test_polys_.encoded.data = data;
        }
        static constexpr size_t n_ = 63;
        static constexpr size_t k_ = 51;
        static constexpr size_t t_ = (n_ - k_) / 6;
        static std::vector <std::bitset <k_>> vector_of_message_polynomials;
        polynomialData<n_, k_> test_polys_{};
};

class Bch6345 {
    public:
        explicit Bch6345(const std::bitset <45> &data) {
            test_polys_.encoded.data = data;
        }
        static constexpr size_t n_ = 63;
        static constexpr size_t k_ = 45;
        static constexpr size_t t_ = (n_ - k_) / 6;
        static std::vector <std::bitset <k_>> vector_of_message_polynomials;
        polynomialData<n_, k_> test_polys_{};
};

class Bch4836 {
    public:
        explicit Bch4836(const std::bitset <36> &data) {
            test_polys_.encoded.data = data;
        }
        static constexpr size_t n_ = 48;
        static constexpr size_t k_ = 36;
        static constexpr size_t t_ = (n_ - k_) / 6;
        static std::vector <std::bitset <k_>> vector_of_message_polynomials;
        polynomialData<n_, k_> test_polys_{};
};

class Bch4830{
    public:
        explicit Bch4830(const std::bitset <30> &data) {
            test_polys_.encoded.data = data;
        }
        static constexpr size_t n_ = 48;
        static constexpr size_t k_ = 30;
        static constexpr size_t t_ = (n_ - k_) / 6;
        static std::vector <std::bitset <k_>> vector_of_message_polynomials;
        polynomialData<n_, k_> test_polys_{};
};

using bchType = std::variant<std::unique_ptr<Bch6351>, std::unique_ptr<Bch6345>, 
                             std::unique_ptr<Bch4836>, std::unique_ptr<Bch4830>>;

namespace bch {
    codeType code_type;

    std::string filename;
    std::mutex Mutex;

    int error_probability;
    template <size_t N>
    mathHelper<N>* bch_math;
    // the big dawg
    std::vector <bchType> BCH_objects;
    std::vector <unsigned char> recovered_charstream;
    std::vector <unsigned char> modified_charstream;
    template <size_t X>
    int MSB(const std::bitset <X> &Polynomial);
    template <size_t X>
    void verbosePolynomial(const std::bitset <X> &Polynomial);
    int multiplyIntPolynomials(int mulitplicand, int multiplicator);
    template <size_t X>
    void reverseBitset(std::bitset <X> &Polynomial, int Shift);
    template<size_t X>
    std::bitset <X> tempReverseBitset(std::bitset <X> Polynomial, int Shift);
    template <size_t X>
    std::pair<std::bitset <X>, std::bitset <X>> divideBitsetPolynomials(
		const std::bitset <X> &dividend, 
		const std::bitset <X> &divisor);
};
