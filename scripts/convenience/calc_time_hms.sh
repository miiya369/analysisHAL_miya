#!/bin/sh

if test $# -ne 3
then
   echo " usage : sh $0 [tiem1] [+/-/x] [time2]"
   echo "         The [time] mast be form with [hour:minutes:second]"
   exit -1
fi

time1_hour=`echo $1 | sed -e "s/^\(.*\):\(.*\):\(.*\)/\1/"`
time1_min=`echo $1 | sed -e "s/^\(.*\):\(.*\):\(.*\)/\2/"`
time1_sec=`echo $1 | sed -e "s/^\(.*\):\(.*\):\(.*\)/\3/"`
time1=$((time1_hour * 60*60 + time1_min * 60 + time1_sec))

if [ $2 = "x" ]
then
    time2=$3
    time=$((time1 * time2))
else
    time2_hour=`echo $3 | sed -e "s/^\(.*\):\(.*\):\(.*\)/\1/"`
    time2_min=`echo $3 | sed -e "s/^\(.*\):\(.*\):\(.*\)/\2/"`
    time2_sec=`echo $3 | sed -e "s/^\(.*\):\(.*\):\(.*\)/\3/"`
    time2=$((time2_hour * 60*60 + time2_min * 60 + time2_sec))
    
    if test $2 = "+"
    then
	time=$((time1 + time2))
    fi
    
    if test $2 = "-"
    then
	time=$((time1 - time2))
    fi    
fi

time_sec=$(((time % (60*60)) % 60))
time_min=$((((time-time_sec)/60) % 60))
time_hour=$(((time - time_min * 60 - time_sec) / (60*60)))

echo "`printf %02d $time_hour`:`printf %02d $time_min`:`printf %02d $time_sec`"
