#include <sys/mman.h>

#include "bch_simulator.hpp"
#include "bch_utils.hpp"
#include "bch_logger.hpp"

#define CURRENT_TIME std::chrono::high_resolution_clock::now();
#define GET_DURATION(a) std::chrono::duration_cast<std::chrono::duration<float>>(a);

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
	const size_t file_byte_size = static_cast<unsigned>(sb.st_size);
	
	char* buffer = static_cast<char*>(mmap(
										NULL, 
										file_byte_size, 
										PROT_READ, MAP_PRIVATE, 
										fd, 
										0));

	std::string file_suffix = bchFirstInit();

	std::string image_with_errors_path = bch::filename.substr(0, bch::filename.length()-4).append("_with_errors_" + file_suffix);
	std::string image_fixed_path = bch::filename.substr(0, bch::filename.length()-4).append("_fixed_" + file_suffix);

	remove(image_with_errors_path.c_str());
	remove(image_fixed_path.c_str());

	std::ofstream image_with_errors (image_with_errors_path.c_str(), std::ios::out | std::ios::app | std::ios::binary);
	std::ofstream image_fixed (image_fixed_path.c_str(), std::ios::out | std::ios::app | std::ios::binary);
	
    // resize vectors:
    const size_t number_of_message_polynomials = resizeMainVectors(file_byte_size);

	std::vector<std::thread> threads;
	const size_t num_threads = std::thread::hardware_concurrency();
	std::cout << "number of detected threads: " << num_threads << std::endl;

	// message_polynomials_thread_group --> count of group of message polynomials that will be used by a single thread
	const size_t message_polynomials_thread_group = (number_of_message_polynomials%num_threads)? 
													number_of_message_polynomials/num_threads :
													number_of_message_polynomials/num_threads - 1;

	// message_bytes_thread_group --> count of group of message bytes that will be used by a single thread
	const size_t message_bytes_thread_group = file_byte_size/num_threads;


	std::cout << "Parsing image file..." << std::endl;
	auto start = CURRENT_TIME;

    divideImageBytesToBitsets(buffer, message_bytes_thread_group, message_polynomials_thread_group,
                                num_threads, file_byte_size, threads);

	auto stop = CURRENT_TIME;
	auto duration = GET_DURATION(stop - start);
    
    auto print_duration = [](auto dur)->std::string { return std::to_string(dur.count()).substr(0, std::to_string(dur.count()).length()-3) + " seconds"; };
    
    std::cout << "Parsing image file done and it took: " << print_duration(duration) << std::endl;

	start = CURRENT_TIME;

    mathStructInit();

	std::cout<<"Please be patient, starting coding and decoding process...\n";

    startMainProcess(number_of_message_polynomials, message_polynomials_thread_group,
                        num_threads, threads);

	stop = CURRENT_TIME;

	auto main_duration = GET_DURATION(stop - start);

	std::cout << "Coding and decoding process done and it took: " << print_duration(main_duration) << std::endl;

	std::cout << "Converting modified and recovered data from bitsets to bytes..." << std::endl;

	start = CURRENT_TIME;

    divideImageBitsetsToBytes(message_bytes_thread_group, message_polynomials_thread_group,
                                num_threads, file_byte_size, threads);

	stop = CURRENT_TIME;
	duration = GET_DURATION(stop - start);

	std::cout << "Converting modified and recovered data from bitsets to bytes done and it took: " << print_duration(duration) << std::endl;
	
	// write bytes to new files and get diff
	image_with_errors.write(buffer, RESERVED_BYTES);
	image_with_errors.write(bch::received_charstream.data() + RESERVED_BYTES, static_cast<signed>(bch::received_charstream.size()) - RESERVED_BYTES);
	image_fixed.write(buffer, RESERVED_BYTES);
	image_fixed.write(bch::decoded_charstream.data() + RESERVED_BYTES, static_cast<signed>(bch::decoded_charstream.size()) - RESERVED_BYTES);

	finalLogsAndCleanup(number_of_message_polynomials, image_with_errors_path, image_fixed_path, main_duration);

	munmap(buffer, file_byte_size);
	close(fd);
	image_with_errors.close();
	image_fixed.close();

	return 0;
}
