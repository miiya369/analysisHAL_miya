#!/bin/sh

if test $# -eq 0
then
   echo " usage : sh $0 [job.out1] [job.out2] ..."
   exit 1
fi

OFILE=($@)
echo " #.job.out file = ${#OFILE[*]}"

count=0
err=0
for i in `seq 0 1 $((${#OFILE[*]}-1))`
do
   cat ${OFILE[$i]} | grep "FATAL ERROR" -A 4
   count=$((count+`grep "GAUGE FIXING END" ${OFILE[$i]} | wc -l`))
   err=$((err+`grep "FATAL ERROR" ${OFILE[$i]} | wc -l`))
done

echo " #.end = $count + $err = $((count+err))"