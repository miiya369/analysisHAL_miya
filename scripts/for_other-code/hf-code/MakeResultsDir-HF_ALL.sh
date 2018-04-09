#!/bin/sh

if [ $# -ne 2 ]; then
    echo " usage : sh `basename $0` [working dir] [conf list]"
    exit 1
fi

working_dir=$1
conf_list=$2

for CONF in `cat $conf_list`; do
    for ch in 0 4; do
	mkdir -p "$working_dir/BBwave.dir.S$ch.00/$CONF" # BBwave: S0 and S4
    done
    
    for ch in 0 1 2 3 4 5 6 7 8; do
	mkdir -p "$working_dir/BBwave.dir.S1.0$ch/$CONF" # BBwave: S1
	mkdir -p "$working_dir/BBwave.dir.S3.0$ch/$CONF" # BBwave: S3
    done
    
    for ch in 00 01 02 03 04 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 31 32 33 35; do
	mkdir -p "$working_dir/BBwave.dir.S2.$ch/$CONF" # BBwave: S2
    done

    for ch in 0 1 2; do
        mkdir -p "$working_dir/PPwave.dir.S$ch.00/$CONF" # PPwave: S0 and S1 and S2
    done

    mkdir -p "$working_dir/Proj.DDwave.dir.S6.00/spin0z+0.0z+0/$CONF" # DDwave: S6.00
    
    for sp in 1 2; do
        mkdir -p "$working_dir/Proj.DOwave.dir.S3.35/spin${sp}_ave/$CONF" # DOwave: S3.35, spin = 1, 2
	mkdir -p "$working_dir/Proj.DOwave.dir.S5.00/spin${sp}_ave/$CONF" # DOwave: S5.00, spin = 1, 2
    done
    
    for had in PS SS multi; do
	mkdir -p "$working_dir/correlator.$had.dir/$CONF" # Correlators
   done
done
