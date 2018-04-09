#!/bin/bash

if [ $# -ne 3 ]; then
   echo "usage: `basename $0` [h:m:s] [+,-,x,/,%] [h:m:s or N(for x,/,%)]"; exit -1
fi

t1h=`echo $1 | cut -d ":" -f 1 | bc`; t2h=`echo $3 | cut -d ":" -f 1 | bc`
t1m=`echo $1 | cut -d ":" -f 2 | bc`; t2m=`echo $3 | cut -d ":" -f 2 | bc`
t1s=`echo $1 | cut -d ":" -f 3 | bc`; t2s=`echo $3 | cut -d ":" -f 3 | bc`
t1=$((3600*t1h + 60*t1m + t1s)); t2=$((3600*t2h + 60*t2m + t2s))

case $2 in
    '+') t=$((t1 + t2));;
    '-') t=$((t1 - t2));;
    'x') t=$((t1 * $3));;
    '/') t=$((t1 / $3));;
    '%') t=$((t1 % $3));;
     * ) echo "ERROR: Invalid operator '$2'"; exit -1
esac

ts=$((( t % 3600) % 60))
tm=$((((t-ts)/60) % 60))
th=$(((t-tm*60-ts)/3600))

echo "`printf %02d $th`:`printf %02d $tm`:`printf %02d $ts`"
