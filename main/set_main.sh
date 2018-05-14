#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: sh `basename $0` [main.cpp]"
    exit
fi

main_cpp=$1

rm -f main_anaHAL.cpp

ln -s $main_cpp main_anaHAL.cpp
