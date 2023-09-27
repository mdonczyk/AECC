#include "bch_utils.hpp"

// create instances for extern namespace vars
namespace bch {
	int error_probability{};
    codeType code_type{};
    std::string filename{};
    std::mutex Mutex{};
}

int parse_arguments(
		const int argc, 
		char* argv[])
{
	int fd = FAIL;
	int c;

	if (argc != 2 && argc != 7 && argc != 8) {
		std::cout << "Invalid number of arguments"<< std::endl;
		return FAIL;
	}

	while ((c = getopt (argc, argv, "i:c:p:vh")) != FAIL) {
		switch (c) {
			case 'h':
				printHelpMessage(argv[0]);
				exit(0);
			case 'i':
				bch::filename = optarg;
				fd = open(optarg, O_RDONLY);
				if (fd == FAIL) {
					std::cout << "Incorrect Image file name" << std::endl;
					return FAIL;
				}
				break;
			case 'c':
				if (atoi(optarg) < 0 || atoi(optarg) > 3) {
					return FAIL;
				}
				bch::code_type = static_cast<codeType>(atoi(optarg));
				break;
			case 'p':
				bch::error_probability = atoi(optarg);
				if (bch::error_probability < 10 || bch::error_probability > 10000000) {
					std::cout << "Error probability should be between (10 and 10000000)" << std::endl;
					return FAIL;
				}
				break;
			case 'v':
				bch_logger::enable_logging = true;
				break;
			case '?':
					std::cout << optopt << std::endl;
				return FAIL;
			default:
				abort();
		}
	}
	return fd;
}

void printHelpMessage(
		const char *file_name)
{
	std::cout << "Usage:\n"
		<< file_name << " [-h] -i image -p err_prob -c code_type [-v]\n\n"
		<< "Options:\n"
		<< "  -i image	Choose one image from images folder, example: images/image2.bmp.\n"
		<< "		  (image has to be in bmp format without compression)\n"
		<< "  -p err_prob	Give probability between 10 and 10000000. A 1 in err_prob error will occur in \n"
		<< "		  the codeword during the simulated transmission through a noisy medium, example: 1000.\n"
		<< "  -c code_type	Choose bch code type:\n"
		<< "			0 - BCH(63, 51, 5)\n"
		<< "			1 - BCH(63, 45, 7)\n"
		<< "			2 - BCH(48, 36, 5)\n"
		<< "			3 - BCH(48, 30, 7)\n"
		<< "Optional arguments:\n"
		<< "  -h	Show this help message.\n"
		<< "  -v	Enable bch_logger encoding and decoding logs which will print out the whole process to the \n"
		<< "	  terminal, is disabled by default. WARNING! This option causes the threads to run sequentially instead \n"
		<< " 	  of in parallel which combined with printing operations to console causes a severe performance degradation.\n";
}
