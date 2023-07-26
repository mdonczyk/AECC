#!/bin/bash
#For most reliable results use the biggest image available

image_path=$1
probabilities=(20 50 80 100 200 500 1000 2000 5000)
programs=(0 1 2 3)

if [[ $# != 1 ]]; then
    echo "Please pass an image path as argument"
    exit 1
fi

for program in ${programs[@]}; do
    echo $'\n'----------------------------------------- >> BCH_logs.txt
    for probability in ${probabilities[@]}; do
        ./bch_simulator -i $image_path -p $probability -c $program
    done
done
