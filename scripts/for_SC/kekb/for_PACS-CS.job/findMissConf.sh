#!/bin/sh

if test $# -ne 3; then
    echo "usage : sh $0 [Tshift] [Ens] [job Number]"
    exit 1
fi

Tshift=$1
Ens=$2
jobN=$3

TSHIFT=`printf %03d $Tshift`
fname=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/pacs-cs_ens$Ens.RHQ/Output/Tshift$TSHIFT/run.out.Tshift${TSHIFT}_$jobN

MissConf=`tail $fname -n 30000 | grep "@@@ reading" | tail -n 1 | cut -d " " -f 3`

echo $MissConf
echo
tail $fname
