#!/bin/bash

### === Change Here === ###
path_of_anaHAL=$HOME/analysisHAL_miya

options="-g -O3 -Wall -fopenmp"
### =================== ###

if [ $# -ne 1 ]; then
    echo "usage: sh `basename $0` [main.cpp]"; exit
fi

main_file=$1
exec_file=`basename $1 .cpp`.x

g++-5 $options -I$path_of_anaHAL/include -L$path_of_anaHAL/lib \
    $main_file -o $exec_file \
    -lpotential -lfield -lanalysis -lyukawa -lyukawaCORE
