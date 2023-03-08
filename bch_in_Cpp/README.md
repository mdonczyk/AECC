Usage:
./bch6351 [-h] -i <image> -p err_prob [-v]

Options:
  -i <file_name>, --image <file_name>   Choose one image from images folder, example: images/image2.bmp.
  -p err_prob, --probability err_prob   Give probability between (1 and 10000000) that an error will occur in
                                          the codeword during a simulated transmission through a noisy medium, example: 1000.
Optional arguments:
  -h, --help    Show this help message.
  -v, --verbose Enable verbose encoding and decoding logs which will print out the whole process to the
                  terminal, is disabled by default. WARNING! This option causes the threads to run sequenitally instead
                  of in parallel which combined with printing operations to console causes a severe performance degradation.

