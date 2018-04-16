#!/bin/bash

path_of_anaHAL=$HOME/analysisHAL_miya

options="-g -O3 -Wall -fopenmp"

if [ $# -ne 1 ]; then
    echo "usage: sh `basename $0` [main.cpp]"; exit
fi

main_file=$1
main_file_bare=`basename $1 .cpp`

g++ $options -I$path_of_anaHAL/include -I$path_of_anaHAL/lib \
    -lanalysis -lfield -lpotential -lyukawa -lyukawaCORE \
    $main_file -o $main_file_bare
