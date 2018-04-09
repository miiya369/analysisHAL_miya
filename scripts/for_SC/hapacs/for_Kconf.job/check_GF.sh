#!/bin/sh

if test $# -eq 0
then
   echo " usage : sh $0 [job.out1] [job.out2] ..."
   exit 1
fi

OFILE=($@)
echo " #.job.out file    = ${#OFILE[*]}"

precision_miss=0
iter_count=0
iter_count_2=0
iter_min=10000000
iter_max=0
count_line=0
for i in `seq 0 1 $((${#OFILE[*]}-1))`
do
   ITER=(`grep "@ converge at" ${OFILE[$i]} | cut -d = -f 2 | cut -d , -f 1`)
   PREC=(`grep "@ converge at" ${OFILE[$i]} | sed -e 's/\(.*\)[-,+]\(.*\)/\2/'`)
   Nline=`grep "@ converge at" ${OFILE[$i]} | wc -l`
   count_line=$((count_line+Nline))

#   echo ${ITER[*]}
#   echo ${PREC[*]}
#   echo "$Nline ${#ITER[*]} ${#PREC[*]}"
#   exit 1

   for l in `seq 0 1 $(($Nline-1))`
   do
      if test ${PREC[$l]} -ne '15'
      then
         precision_miss=$(( precision_miss + 1 ))
      fi      
      if test ${ITER[$l]} -lt $iter_min
      then
         iter_min=${ITER[$l]}
      fi
      if test ${ITER[$l]} -gt $iter_max
      then
	 iter_max=${ITER[$l]}
      fi
      iter_count=$(( iter_count + ${ITER[$l]} ))
      iter_count_2=$(( iter_count_2 + ${ITER[$l]}*${ITER[$l]} ))
   done
done

echo " #.Miss precision  = $precision_miss"
echo " Minimum iteration = $iter_min"
echo " Maximum iteration = $iter_max"
iter_ave=$(( iter_count / ($Nline*${#OFILE[*]}) ))
iter_ave_2=$(( iter_count_2 / ($Nline*${#OFILE[*]}) ))
iter_err=`echo "sqrt(($iter_ave_2 - $iter_ave*$iter_ave) / ($Nline*${#OFILE[*]} - 1))" | bc`
echo " Average iteration = $iter_ave +/- $iter_err"
echo " #.line = $count_line (par #.job.out = `echo $count_line/${#OFILE[*]} | bc -l`)"