#!/bin/sh

if [ $# -eq 0 ]; then
    echo "usage: sh `basename $0` [src t] ..."
    exit 1
fi

Base=/gwfefs/data/G16030/miyamoto/work.K/work.run-i.traj_mod10_5/run.t_mod24_??/run.t_???.rot_tot/logs.t_???
run_out_lst=/gwfefs/home/tkymiya/scripts/run.out.lst
Check=/gwfefs/home/tkymiya/scripts/CheckJob.CalcConf.py

ConfNs=(0 0 0 0 0 0 0)
OKconf=(0 0 0 0 0 0 0)

for srcT in $@; do
    SRCT=`printf %03d $srcT`
    if [ $? -ne 0 ]; then
	exit 1
    fi
    for ofile in `cat $run_out_lst | sed "s/SRCT/$SRCT/g"`; do
	count=0
	i=0
	Sconf=`echo $ofile | cut -d "." -f 5 | cut -d "-" -f 1`
	Econf=`echo $ofile | cut -d "." -f 5 | cut -d "-" -f 2`
	for conf in `seq $Sconf 10 $Econf`; do
	    ConfNs[$i]=$conf
	    OKconf[$i]=0
	    i=$((i+1))
	done
	for res in `$Check $Base/$ofile`; do
	    count=$((count+1))
	    
	    if [ $count -eq 1 ]; then
                srcT=$res
            elif [ $count -eq 2 ]; then
                conf=$res
            elif [ $count -eq 3 ]; then
		txyz=$res
	    elif [ $count -eq 4 ]; then
		count=0
		for j in `seq 0 1 6`; do
		    if [ ${ConfNs[$j]} -eq $conf -a $res == "Successful" ]; then
			OKconf[$j]=1
		    elif [ ${ConfNs[$j]} -eq $conf -a $res == "Failed" ]; then
			OKconf[$j]=0
		    fi
		done
	    fi
	done
	for j in `seq 0 1 6`; do
            if [ ${OKconf[$j]} -eq 0 ]; then
                echo "$srcT ${ConfNs[$j]} $txyz"
            fi
        done
    done
done
