#!/bin/sh

### Auther: Takaya Miyamoto 
### Date  : Fri Sep  2 00:51:57 JST 2016
###
### Brief : Make the text file for gnuplot
###       : It can be used for Not-Compressed Spin-Projected NBS data.
###

if [ $# -ne 1 ]; then
    echo "usage: `basename $0` [NBS file]"
    exit 1
fi

NBSfile=$1
CrdSize=16
ShowBinary=$HOME/analysis_code-set/scripts/convenience/ShowBinary_Double.py

for iz in `seq 0 1 $((CrdSize-1))`; do
    for iy in `seq 0 1 $((CrdSize-1))`; do
	for ix in `seq 0 1 $((CrdSize-1))`; do
	    if [ $ix -le $((CrdSize/2)) -a $iy -le $ix -a $iz -le $iy ]; then
		index0=$((0+2*(ix+$CrdSize*(iy+$CrdSize*iz))))
		index1=$((1+2*(ix+$CrdSize*(iy+$CrdSize*iz))))
		ir=`echo "sqrt($ix*$ix+$iy*$iy+$iz*$iz)" | bc -l`
		res0=`$ShowBinary $NBSfile 1 $index0 | awk '{print $3}'`
		res1=`$ShowBinary $NBSfile 1 $index1 | awk '{print $3}'`
		echo "$ir $res0 $res1"
	    fi
	done
    done
done
