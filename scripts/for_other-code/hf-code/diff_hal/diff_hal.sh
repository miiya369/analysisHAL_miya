#!/bin/sh

if [ $# -ne 5 ]; then
   echo "usage: `basename $0` [dir list] [conf list] [file list (For corr)] [results1] [results2]"
   exit 1
fi

######### Change Here #########
diff_hal=$HOME/analysis_code-set/scripts/for_other-code/hf-code/diff_hal/diff_hal.py
criterion_err=1e-5
criterion_rate=0.3

t_min=6  # For NBS
t_max=18 # For NBS
###############################

dir_lst=$1
conf_lst=$2
file_lst=$3
results1=$4
results2=$5

echo " ========= Setup ========="
echo " current   dir  = `pwd`"
echo " channel   list = $dir_lst"
echo " conf      list = $conf_lst"
echo " corr file list = $file_lst"
echo " results 1      = $results1"
echo " results 2      = $results2"
echo
echo " The criterion for the order of relative error = $criterion_err"
echo " The criterion for the rate  of #.components   = $criterion_rate"
echo
echo " Time slice = $t_min - $t_max"
echo " ========================="
echo
echo "JOB Start: `date`"
for Dir in `cat $dir_lst`; do
   Flag=""

   if [ `echo $Dir | grep BB | wc -l` -eq 1 ]; then
      Ftype=0
      Fbase=NBSwave.+REP_TIME+???.000.000.000.REP_CONF.NUC_CG05.NUC_CG05
      case $Dir in
         "BBwave.dir.S2.00"|"BBwave.dir.S2.01"|"BBwave.dir.S2.02"|"BBwave.dir.S2.03"|"BBwave.dir.S2.04"|\
         "BBwave.dir.S2.06"|"BBwave.dir.S2.10"|"BBwave.dir.S2.12"|"BBwave.dir.S2.16"|"BBwave.dir.S2.18"|\
         "BBwave.dir.S2.22"|"BBwave.dir.S2.24"|"BBwave.dir.S2.25"|"BBwave.dir.S2.26"|"BBwave.dir.S2.27"|\
         "BBwave.dir.S2.28") Flag="-spin0";;
         "BBwave.dir.S2.23"|"BBwave.dir.S2.33") Flag="-spin1";;
      esac
   elif [ `echo $Dir | grep PP | wc -l` -eq 1 ]; then
      Ftype=1
      Fbase=NBSwave.+REP_TIME+???.000.000.000.REP_CONF.PP_CG05.PP_CG05
   elif [ `echo $Dir | grep DD | wc -l` -eq 1 ]; then
      Ftype=1
      Fbase=NBSwave.+REP_TIME+???.000.000.000.REP_CONF.DD_NR.DD_NR
   elif [ `echo $Dir | grep DO | wc -l` -eq 1 ]; then
      Ftype=1
      Fbase=NBSwave.+REP_TIME+???.000.000.000.REP_CONF.DO_NR.DO_NR
   elif [ `echo $Dir | grep corr | wc -l` -eq 1 ]; then
      Ftype=2
   else
      echo "Invalid directory name; $Dir."
      exit 1
   fi

   for Conf in `cat $conf_lst`; do

      if [ $Ftype != 2 ]; then
         for time in `seq $t_min 1 $t_max`; do
            Time=`printf %03d $time`
            File=`echo $Fbase | sed -e s/REP_TIME/$Time/g`
            File=`echo $File  | sed -e s/REP_CONF/$Conf/g`

            Path1=$results1/$Dir/$Conf/$File
            Path2=$results2/$Dir/$Conf/$File

            for Ndiff in `$diff_hal $Path1 $Path2 $Ftype -diff_rate $criterion_err $Flag`; do
               if [ `awk "BEGIN {if ($Ndiff > $criterion_rate) print 1; else print 0}"` == 1 ]; then
                  echo "$Dir/$Conf/$File"
                  echo "in the $results1 and $results2 are differ; Rate = $Ndiff."
               fi
            done
         done
      else
         for File in `cat $file_lst`; do
            Path1=$results1/$Dir/$Conf/$File
            Path2=$results2/$Dir/$Conf/$File

            for Ndiff in `$diff_hal $Path1 $Path2 $Ftype -diff_rate $criterion_err $Flag`; do
               if [ `awk "BEGIN {if ($Ndiff > $criterion_rate) print 1; else print 0}"` == 1 ]; then
                  echo "$Dir/$Conf/$File"
                  echo "in the $results1 and $results2 are differ; Rate = $Ndiff."
               fi
            done
         done
      fi
   done
   echo "$Dir End : `date`"
done
echo "JOB  End : `date`"
