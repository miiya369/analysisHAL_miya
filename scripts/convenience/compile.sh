#!/bin/sh

if [ $# -ne 1 ]; then
    echo "usage: `basename $0` [C/C++ file]"
    exit 1
fi

CXX=g++
OPTION="-Wall -O3 -g"

codeName=$1
binName=${codeName%.*}

$CXX $OPTION -o $binName $codeName
