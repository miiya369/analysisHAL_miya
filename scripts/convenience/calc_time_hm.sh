#!/bin/sh

if test $# -ne 3
then
   echo " usage : sh $0 [tiem1] [+/-] [time2]"
   echo "         The [time] mast be form with [hour:minutes]"
   exit -1
fi

time1_hour=`echo $1 | sed -e "s/^\(.*\):\(.*\)/\1/"`
time2_hour=`echo $3 | sed -e "s/^\(.*\):\(.*\)/\1/"`
time1_min=`echo $1 | sed -e "s/^\(.*\):\(.*\)/\2/"`
time2_min=`echo $3 | sed -e "s/^\(.*\):\(.*\)/\2/"`

time1=$((time1_hour * 60 + time1_min))
time2=$((time2_hour * 60 + time2_min))

if test $2 = "+"
then
   time=$((time1 + time2))
fi

if test $2 = "-"
then
   time=$((time1 - time2))
fi

time_min=$((time % 60))
time_hour=$(((time - time_min) / 60))

echo "`printf %02d $time_hour`:`printf %02d $time_min`"
