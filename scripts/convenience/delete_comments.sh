#!/bin/sh

# Author: Takaya Miyamoto
# Date  : Fri Sep 16 13:18:08 JST 2016
# Brief : Delete the line contains certain words in some files.

if [ $# -le 1 ]; then
    echo "usage: `basename $0` [Word(s)] [file1] [file2] ..."
    exit 1
fi

args=($@)
Narg=$#
swc=0

echo "Delete the line contains '${args[0]}'. OK? :"; read yn

if [ "a${yn}a" == "aa" ]; then exit 1; fi   # Safety Device
if [ $yn != "y" -a $yn != "yes" -a $yn != "yeah" -a $yn != "yap" ]; then exit 1; fi

for i in `seq 1 1 $((Narg-1))`; do
    sed -i.bak -e "/${args[0]}/d" ${args[$i]}
done
