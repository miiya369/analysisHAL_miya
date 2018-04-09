#!/bin/sh

if [ $# -eq 0 ]; then
    echo "usage: sh `basename $0` [ofile1] [ofile2] ..."
    exit 1
fi

mean=0
mean2=0

for ofile in $@; do    
    Stime=`head -n 1 $ofile | cut -d " " -f 4`
    Etime=`tail -n 1 $ofile | cut -d " " -f 4`
    
    Stime_h=`echo $Stime | cut -d ":" -f 1`
    Stime_m=`echo $Stime | cut -d ":" -f 2`
    Stime_s=`echo $Stime | cut -d ":" -f 3`
    Etime_h=`echo $Etime | cut -d ":" -f 1`
    Etime_m=`echo $Etime | cut -d ":" -f 2`
    Etime_s=`echo $Etime | cut -d ":" -f 3`
    
    Stime_sec=$((10#$Stime_h * 60*60 + 10#$Stime_m * 60 + 10#$Stime_s))
    Etime_sec=$((10#$Etime_h * 60*60 + 10#$Etime_m * 60 + 10#$Etime_s))
    
    if [ $Stime_sec -lt $Etime_sec ]; then
	Elapse_sec=$((Etime_sec - Stime_sec))
    else
	tmp_sec=$((86400 - Stime_sec))
	Elapse_sec=$((tmp_sec + Etime_sec))
    fi
    
    mean=$((mean + Elapse_sec))
    mean2=$((mean2 + Elapse_sec*Elapse_sec))
done

mean=$((mean / $#))
mean2=$((mean2 / $#))

if [ $# -eq 1 ]; then
    err=0
else
    err=`echo "sqrt(($mean2 - $mean*$mean) / ($# - 1))" | bc`
fi

mean_s=$(((mean % 3600) % 60))
mean_m=$((((mean-mean_s)/ 60) % 60))
mean_h=$(((mean-mean_m*60-mean_s) / 3600))

err_s=$(((err % 3600) % 60))
err_m=$((((err-err_s)/ 60) % 60))
err_h=$(((err-err_m*60-err_s) / 3600))

echo "$mean_h:$mean_m:$mean_s +/- $err_h:$err_m:$err_s ($#)"
