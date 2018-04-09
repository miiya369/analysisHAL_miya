#!/bin/sh

cd $HOME/report/performance

f=`ls *.end`
highest=0.0

for F in $f
do
    tmp=(`grep 'GFLOPS(MA) per node:' $F`)
    TMP=${tmp[3]}

    echo $TMP | grep "e" > /dev/null 2>&1

    if [ $? -ne 0 ]
    then
	if [ `echo "$TMP > $highest" | bc` -eq 1 ]
	then
            highest=$TMP
            Fname=$F
	fi
    fi
done

echo "HIGHEST JOB: $Fname"
echo "GFLOPS/NODE: $highest"
