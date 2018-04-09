#!/bin/sh

if [ $# -ne 3 ]; then
   echo "usage: `basename $0` [Job ID] [Ens] [Job queue]"
   exit 1
fi

JobID=$1
Ens=$2
Que=$3
Que_org=`bllq -l $JobID | grep 'Class:' | awk '{print $2}'`
JobFile=`bllq -l $JobID | grep 'Cmd:' | awk '{print $2}' | cut -d "/" -f 7 | cut -d "." -f 2-3`
Tshift=`echo $JobFile | cut -d "_" -f 1`

File=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/pacs-cs_ens$Ens.RHQ/Jobscript/$Tshift/$JobFile

if [ ! -e $File ]; then
    echo "ERROR: No such a file, $JobFile"
    exit 1
fi

echo "Change $File"
echo "$Que_org -> $Que"

cat $File | sed s/$Que_org/$Que/g > ${File}_chq
bllcancel $JobID
bllsubmit ${File}_chq
