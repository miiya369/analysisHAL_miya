#!/bin/sh

if test $# -lt 9
then
   echo " usage : sh $0 [job.out1] [job.out2] ..."
   echo " Need to #.job.out >= 9"
   exit 1
fi

OFILE=($@)
echo " #.job.out file    = ${#OFILE[*]}"

iter_count=0
iter_count_2=0
iter_min=10000000
iter_max=0

width_hist=50

dev_hist=$((1000000/width_hist))
high_hist=$((${#OFILE[*]} / 9))
for i in `seq 0 1 $((width_hist-1))`
do
   iter_hist+=("0")
done

for i in `seq 0 1 $((${#OFILE[*]}-1))`
do
   ITER=(`grep "@ converge at" ${OFILE[$i]} | cut -d = -f 2 | cut -d , -f 1`)
   Nline=`grep "@ converge at" ${OFILE[$i]} | wc -l`

#   echo ${ITER[*]}
#   echo ${PREC[*]}
#   echo "$Nline ${#ITER[*]} ${#PREC[*]}"
#   exit 1

   for l in `seq 0 1 $(($Nline-1))`
   do
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
      iter_hist[$((${ITER[$l]}/$dev_hist))]=$((iter_hist[$((${ITER[$l]}/$dev_hist))] + 1))
   done
done

echo " Minimum iteration = $iter_min"
echo " Maximum iteration = $iter_max"
iter_ave=$(( iter_count / ($Nline*${#OFILE[*]}) ))
iter_ave_2=$(( iter_count_2 / ($Nline*${#OFILE[*]}) ))
iter_err=`echo "sqrt(($iter_ave_2 - $iter_ave*$iter_ave) / ($Nline*${#OFILE[*]} - 1))" | bc`
echo " Average iteration = $iter_ave +/- $iter_err"

echo
echo " Histogram for #.iter (\"|\" = $high_hist)"
for i in `seq 0 1 $((width_hist-1))`
do
   echo -n " `printf %06d $((dev_hist * $i))` : "
   if test ${iter_hist[$i]} -lt $high_hist
   then
      echo -n "${iter_hist[$i]}"
   else
      for j in `seq 1 $high_hist ${iter_hist[$i]}`
      do
         echo -n "|"
      done
   fi
   echo
done
