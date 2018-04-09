#!/bin/sh

if test $# -ne 3
then 
    echo " usage : sh $0 [depend job num 1] [depend job num 2] [job script]"
    exit -1
fi

qsub -W depend=afterok:$1:$2 $3
