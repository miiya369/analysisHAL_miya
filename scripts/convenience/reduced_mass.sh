#!/bin/sh

if [ $# -ne 2 ]; then
   echo " usage : sh `basename $0` [mass1] [mass2]"
   exit 1
fi

mass1=$1
mass2=$2

echo "scale=6; ($mass1 * $mass2) / ($mass1 + $mass2)" | bc
