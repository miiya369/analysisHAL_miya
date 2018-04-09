#!/bin/sh

if [ $# -ne 3 ]
then
    echo "usage: sh $0 [Path] [Min file num] [Max file num]"
    exit 1
fi

path=$1
min=$2
max=$3
files=(`ls $1`)

for i in `seq $min 1 $max`
do
    echo "$path/${files[$i-1]}"
done
