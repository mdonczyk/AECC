#include <iostream>
#include <memory>
#include <iomanip>

#include "bch_utils.hpp"

// create instances for extern namespace vars
namespace bch {
	int error_probability{};
    codeType code_type{};
    std::string filename{};
    std::mutex Mutex{};
}

void printHelpMessage(
		const char *file_name)
{
	std::cout << "Usage:\n"
		<< file_name << " [-h] -i image -p err_prob -c code_type [-v]\n\n"
		<< "Options:\n"
		<< "  -i image			Choose one image from images folder, example: images/image2.bmp.\n"
		<< "  -p err_prob		Give probability between (10 and 10000000) that a 1 in err_prob error will occur in \n"
		<< "			   			the codeword during a simulated transmission through a noisy medium, example: 1000.\n"
		<< "  -c code_type		Choose bch code type:\n"
		<< "						0 - BCH(63, 51, 5)\n"
		<< "						1 - BCH(63, 45, 7)\n"
		<< "						2 - BCH(48, 36, 5)\n"
		<< "						3 - BCH(48, 30, 7)\n"
		<< "Optional arguments:\n"
		<< "  -h	Show this help message.\n"
		<< "  -v	Enable verbose_flag encoding and decoding logs which will print out the whole process to the \n"
		<< "	   terminal, is disabled by default. WARNING! This option causes the threads to run sequentially instead \n"
		<< " 	   of in parallel which combined with printing operations to console causes a severe performance degradation.\n";
}
