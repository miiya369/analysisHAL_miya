#!/bin/sh

if [ $# -lt 5 ]; then
   echo "usage: `basename $0` [working dir] [base name of confs] [#.org conf] [bin size] [dir name 1] [dir name 2] ..."
   exit 1
fi

working_dir=$1
directories=($@)
nbase=$2
Nconf=$3
bin_size=$4

if [ $((Nconf%bin_size)) -ne 0 ]; then
   echo "ERROR: #.conf cannot be divided by bin size."
   exit 1
fi

N_bin=$((Nconf/bin_size))

for ibin in `seq 0 1 $((N_bin-1))`; do
    CONF=${nbase}_bin`printf %03d $bin_size`-`printf %03d $ibin`
    
    for i in `seq 5 1 $#`; do
	if [ ! -e $working_dir/${directories[i-1]}/$CONF ]; then
            mkdir -p $working_dir/${directories[i-1]}/$CONF
	fi
    done
done
