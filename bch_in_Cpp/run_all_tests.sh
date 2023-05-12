#!/bin/bash
#For most reliable results use the biggest image available

image_path=$1

echo $'\n'----------------------------------------- >> BCH_logs.txt

./bch6351 -i $image_path -p 20
./bch6351 -i $image_path -p 50
./bch6351 -i $image_path -p 80
./bch6351 -i $image_path -p 100
./bch6351 -i $image_path -p 200
./bch6351 -i $image_path -p 500
./bch6351 -i $image_path -p 1000
./bch6351 -i $image_path -p 2000
./bch6351 -i $image_path -p 5000

echo $'\n'----------------------------------------- >> BCH_logs.txt

./bch6345 -i $image_path -p 20
./bch6345 -i $image_path -p 50
./bch6345 -i $image_path -p 80
./bch6345 -i $image_path -p 100
./bch6345 -i $image_path -p 200
./bch6345 -i $image_path -p 500
./bch6345 -i $image_path -p 1000
./bch6345 -i $image_path -p 2000
./bch6345 -i $image_path -p 5000

echo $'\n'----------------------------------------- >> BCH_logs.txt

./bch4836 -i $image_path -p 20
./bch4836 -i $image_path -p 50
./bch4836 -i $image_path -p 80
./bch4836 -i $image_path -p 100
./bch4836 -i $image_path -p 200
./bch4836 -i $image_path -p 500
./bch4836 -i $image_path -p 1000
./bch4836 -i $image_path -p 2000
./bch4836 -i $image_path -p 5000

echo $'\n'----------------------------------------- >> BCH_logs.txt

./bch4830 -i $image_path -p 20
./bch4830 -i $image_path -p 50
./bch4830 -i $image_path -p 80
./bch4830 -i $image_path -p 100
./bch4830 -i $image_path -p 200
./bch4830 -i $image_path -p 500
./bch4830 -i $image_path -p 1000
./bch4830 -i $image_path -p 2000
./bch4830 -i $image_path -p 5000

