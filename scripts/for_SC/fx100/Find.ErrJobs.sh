#!/bin/sh

if [ $# -eq 0 ]; then
    echo "usage: sh `basename $0` [src t] ..."
    exit 1
fi

Base=/gwfefs/data/G16030/miyamoto/work.K/work.run-i.traj_mod10_5/run.t_mod24_??/run.t_???.rot_tot/logs.t_???
run_out_lst=/gwfefs/home/tkymiya/scripts/run.out.lst

for srcT in $@; do
    SRCT=`printf %03d $srcT`
    if [ $? -ne 0 ]; then
	exit 1
    fi
    for ofile in `cat $run_out_lst | sed "s/SRCT/$SRCT/g"`; do
	if [ ! -e $Base/$ofile ]; then
	    echo "$ofile : Not exist"
	    continue
	fi
	
	tmp0=`grep ERROR $Base/$ofile | wc -l`
	if [ $tmp0 -ne 0 ]; then
            echo "$ofile : Error happened"
            continue
	fi

	tmp4=`ls -l $Base/$ofile | cut -d " " -f 5`
        if [ $tmp4 -lt 100000 ]; then
            echo "$ofile : Invalid file size : Too small file size"
            continue
        fi
	
	tmp1=`tail -n 1 $Base/$ofile | cut -d " " -f 4`
	tmp2=`echo $tmp1 | cut -d ":" -f 1`
	expr $tmp2 + 1 > /dev/null 2>&1
	if [ $? -ne 0 ]; then
	    echo "$ofile : Time limit exceeded"
	    continue
	fi
	
	tmp3=`ls -l $Base/$ofile | cut -d " " -f 5`
	if [ $tmp3 -lt 705000 -o $tmp3 -gt 715000 ]; then
	    echo -n "$ofile : Invalid file size"
	    if [ `grep "restarts" $Base/$ofile | wc -l` -ne 0 ]; then
		echo " : restarted."
	    else
		echo " : Something happened"
	    fi
	    continue
	fi
    done
done
