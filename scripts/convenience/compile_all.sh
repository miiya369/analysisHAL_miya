#!/bin/sh

EXT=cpp
CXX=g++
OPTION="-Wall -O3 -g"

for code in *.$EXT; do
    codeName=$code
    binName=${codeName%.*}
    echo "Compile: $codeName ..."
    $CXX $OPTION -o $binName $codeName
done
