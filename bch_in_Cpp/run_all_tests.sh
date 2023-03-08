#!/bin/bash

err_prob=$1

./bch6351 -i images/lenna.bmp -p $err_prob
./bch6351 -i images/image1.bmp -p $err_prob
./bch6351 -i images/image2.bmp -p $err_prob
./bch6351 -i images/image3.bmp -p $err_prob
./bch6351 -i images/image4.bmp -p $err_prob

./bch6345 -i images/lenna.bmp -p $err_prob
./bch6345 -i images/image1.bmp -p $err_prob
./bch6345 -i images/image2.bmp -p $err_prob
./bch6345 -i images/image3.bmp -p $err_prob
./bch6345 -i images/image4.bmp -p $err_prob

./bch4836 -i images/lenna.bmp -p $err_prob
./bch4836 -i images/image1.bmp -p $err_prob
./bch4836 -i images/image2.bmp -p $err_prob
./bch4836 -i images/image3.bmp -p $err_prob
./bch4836 -i images/image4.bmp -p $err_prob

./bch4830 -i images/lenna.bmp -p $err_prob
./bch4830 -i images/image1.bmp -p $err_prob
./bch4830 -i images/image2.bmp -p $err_prob
./bch4830 -i images/image3.bmp -p $err_prob
./bch4830 -i images/image4.bmp -p $err_prob
