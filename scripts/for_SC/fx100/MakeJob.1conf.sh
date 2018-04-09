#!/bin/sh

if [ $# -ne 3 ]; then
    echo "usage: sh `basename $0` [src T] [conf] [txyz]"
    exit 1
fi

SRCT=$1
conf=$2
txyz=$3

JobScript=/gwfefs/home/tkymiya/scripts/PrintJobScript.1conf.sh
Parameter=/gwfefs/home/tkymiya/scripts/PrintParameter.1conf.sh

Base=/gwfefs/data/G16030/miyamoto/work.K/work.run-i.traj_mod10_5
JobScriptDir=Jobs.dir.run.v04.for_FX100.run-i.doi.2016-04-17/Jobs
ParameterDir=ParamsHadron.dir.ud_e-06.s_e-12.v03.for_FX100.run-i.doi.2015-06-14/Params

if [ $SRCT == "002" -o $SRCT == "026" -o $SRCT == "050" -o $SRCT == "074" ]; then
    Mod24=02
elif [ $SRCT == "006" -o $SRCT == "030" -o $SRCT == "054" -o $SRCT == "078" ]; then
    Mod24=06
else
    echo "ERROR: Invalid src T"
    exit 1
fi

CONF=`printf %06d $conf`

$JobScript $SRCT $conf $txyz > $Base/$JobScriptDir.t_mod24_0$Mod24/job.t_$SRCT.$CONF.rot_$txyz
$Parameter $SRCT $conf $txyz > $Base/$ParameterDir.t_mod24_0$Mod24/ParamsHadron.t_$SRCT.$CONF.rot_$txyz
