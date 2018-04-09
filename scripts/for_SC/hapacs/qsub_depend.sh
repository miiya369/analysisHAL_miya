#!/bin/sh

if test $# -ne 2
then 
    echo " usage : sh $0 [depend job num] [job script]"
    exit -1
fi

qsub -W depend=afterok:$1 $2
