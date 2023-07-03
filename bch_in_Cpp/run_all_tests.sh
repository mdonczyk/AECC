#!/bin/bash
#For most reliable results use the biggest image available

image_path=$1
probabilities=(20 50 80 100 200 500 1000 2000 5000)
programs=(bch6351 bch6345 bch4836 bch4830)

if [[ $# -eq 0 || $# -gt 1 ]]; then
    echo "Please pass an image path as argument"
    exit 1
fi

for program in ${programs[@]}; do
    echo $'\n'----------------------------------------- >> BCH_logs.txt
    for probability in ${probabilities[@]}; do
        ./$program -i $image_path -p $probability
    done
done
