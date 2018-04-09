#!/bin/sh

if [ $# -ne 6 ]
then
    echo "usage: sh $0 [] [] ... [] <- put 6 characters for safety."
    exit 1
fi

Jobscripts=`ls`
for JOB in $Jobscripts
do
    bllsubmit $JOB
done
