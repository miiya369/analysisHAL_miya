#!/bin/sh

if [ $# -ne 3 ]
then
    echo "usage: sh $0 [Path to Gfiles] [#.conf/job] [Job Num.]"
    exit 1
fi

PathToConf=$1
Nconf=`ls $PathToConf | wc -w`
confs=(`ls $1`)

Nconf_job=$2
jobN=$3

confN_s=$((jobN    * Nconf_job + 1))
confN_e=$((confN_s + Nconf_job - 1))

if [ $confN_e -gt $Nconf ]
then
    confN_e=$Nconf
fi

for i in `seq $confN_s 1 $confN_e`
do
    echo "$PathToConf/${confs[$i-1]}"
done
