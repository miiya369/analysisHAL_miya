#!/bin/sh

if [ $# -ne 2 ]; then
    echo "usage: `basename $0` [ifname] [mass]"
    exit 1
fi

###### Change Here ######
#SolveSch=/$HOME/analysis_code-set/scripts/convenience/Solve.Sch.Gauss_exp.py
SolveSch=/Users/miiya/Dropbox/programs/analysis_code-set/scripts/convenience/Solve.Sch.Gauss_exp.py

lat_space=0.0907

max_r=10
Nret=10
Nbase=20
#########################

ifile=$1
mass=$2
tmp_d=(999999.0 0.0)
count=0

while :; do
    EigE=(`$SolveSch $ifile -mass $mass -lat_space $lat_space -max_r $max_r -Nbase $Nbase -stochastic 2>&1 \
	| grep "Eigen energy" \
	| cut -d " " -f 5-7`)
    
    if [ ${#EigE[@]} -eq 0 ]; then continue; fi
    if [ ${#EigE[@]} -ne 3 ]; then
	echo "Something Happend, exit."; exit 1
    fi
    
    if [ `echo "${EigE[0]} < ${tmp_d[0]}" | bc` -eq 1 ]; then
	count=$((count+1))
	tmp_d[0]=${EigE[0]}
	tmp_d[1]=${EigE[2]}
	echo "calculation[`printf %03d $count`]: E = ${tmp_d[0]} +/- ${tmp_d[1]}"
    fi
    
    if [ $count -eq $Nret ]; then break; fi
done
