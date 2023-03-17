#!/bin/bash

image_path=$1

./bch6351 -i $image_path -p 10
./bch6351 -i $image_path -p 50
./bch6351 -i $image_path -p 100
./bch6351 -i $image_path -p 200
./bch6351 -i $image_path -p 500
./bch6351 -i $image_path -p 1000
./bch6351 -i $image_path -p 2000
./bch6351 -i $image_path -p 5000

./bch6345 -i $image_path -p 10
./bch6345 -i $image_path -p 50
./bch6345 -i $image_path -p 100
./bch6345 -i $image_path -p 200
./bch6345 -i $image_path -p 500
./bch6345 -i $image_path -p 1000
./bch6345 -i $image_path -p 2000
./bch6345 -i $image_path -p 5000

./bch4836 -i $image_path -p 10
./bch4836 -i $image_path -p 50
./bch4836 -i $image_path -p 100
./bch4836 -i $image_path -p 200
./bch4836 -i $image_path -p 500
./bch4836 -i $image_path -p 1000
./bch4836 -i $image_path -p 2000
./bch4836 -i $image_path -p 5000

./bch4830 -i $image_path -p 10
./bch4830 -i $image_path -p 50
./bch4830 -i $image_path -p 100
./bch4830 -i $image_path -p 200
./bch4830 -i $image_path -p 500
./bch4830 -i $image_path -p 1000
./bch4830 -i $image_path -p 2000
./bch4830 -i $image_path -p 5000

