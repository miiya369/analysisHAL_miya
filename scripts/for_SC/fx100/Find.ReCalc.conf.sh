#!/bin/sh

if [ $# -eq 0 ]; then
    echo "usage: sh `basename $0` [src t] ..."
    exit 1
fi

Base=/gwfefs/data/G16030/miyamoto/work.K/work.run-i.traj_mod10_5/run.t_mod24_??/run.t_???.rot_tot/recalc.job.logs

for srcT in $@; do
    SRCT=`printf %03d $srcT`
    if [ $? -ne 0 ]; then
	exit 1
    fi
    List=(`ls -l $Base/run.out.run_i.t_$SRCT.*.rot_?`)
    if [ $? -eq 0 ]; then
	for i in `seq 9 10 ${#List[@]}`; do
	    tmp1=`echo ${List[$i]} | cut -d "." -f 12 | cut -d "_" -f 2`
	    tmp2=`echo ${List[$i]} | cut -d "." -f 13`
	    tmp3=`echo ${List[$i]} | cut -d "." -f 14 | cut -d "_" -f 2`
	    Tmp2=`expr $tmp2 + 0`
	    echo "$tmp1 $Tmp2 $tmp3"
	done
    fi
done
