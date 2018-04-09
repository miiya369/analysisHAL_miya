#!/bin/sh

if [ $# -ne 5 ]; then
    echo "usage : sh $0 [Tshift] [Ens] [job Number] [Miss Conf] [Que]"
    exit 1
fi

Tshift=$1
TSHIFT=`printf %03d $Tshift`
Ens=$2
jobN=$3
MissConf=$4
Que=$5

PrintJobscript=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/script/PrintJobscript_Bridge.sh
FindMissConf=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/script/findMissConf.sh
BASE1=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/pacs-cs_ens$Ens.RHQ/Jobscript/Tshift$TSHIFT
BASE2=/bghome/scnfqcd/tkymiya/bridge-1.3.2-hal.v01/src2
BASE3=/bghome/scnfqcd/tkymiya/gfile_list/pacs-cs_ens${Ens}_gfix

wdir=`cat $BASE1/Tshift${TSHIFT}_${jobN}.job | grep workingdir | cut -d "=" -f 2`
JobName=`cat $BASE1/Tshift${TSHIFT}_${jobN}.job | grep jobname | cut -d "=" -f 2`
ofile=`cat $BASE1/Tshift${TSHIFT}_${jobN}.job | grep outputfile | cut -d "=" -f 2`_avoid.$MissConf
efile=`cat $BASE1/Tshift${TSHIFT}_${jobN}.job | grep errorfile | cut -d "=" -f 2`_avoid.$MissConf
param=`cat $BASE1/Tshift${TSHIFT}_${jobN}.job | grep arguments | cut -d " " -f 2`
gfile=`cat $BASE1/Tshift${TSHIFT}_${jobN}.job | grep arguments | cut -d " " -f 4`_avoid.$MissConf

$PrintJobscript $Que $wdir $JobName $ofile $efile $param $gfile \
    > $BASE1/Tshift${TSHIFT}_${jobN}_avoid.$MissConf.job

if [ ! -e $BASE2/run.remain_pacs-cs.ens$Ens.Tshift$TSHIFT ]; then
    cp -r $BASE2/run.remain_pacs-cs.ens$Ens.Tshift000 $BASE2/run.remain_pacs-cs.ens$Ens.Tshift$TSHIFT
fi
cd $BASE2/run.remain_pacs-cs.ens$Ens.Tshift$TSHIFT
LinkToResults.sh $Ens $Tshift
wdir2=`pwd`

ofile2=`cat $BASE1/Tshift${TSHIFT}_${jobN}.job | grep outputfile | cut -d "=" -f 2 | cut -d "_" -f 1-3`_$MissConf
efile2=`cat $BASE1/Tshift${TSHIFT}_${jobN}.job | grep errorfile | cut -d "=" -f 2 | cut -d "_" -f 1-3`_$MissConf
gfile2=`cat $BASE1/Tshift${TSHIFT}_${jobN}.job | grep arguments | cut -d " " -f 4 | cut -d "_" -f 1-7`_$MissConf

$PrintJobscript 32m $wdir2 RemainJob $ofile2 $efile2 $param $gfile2 \
    > $BASE1/Tshift${TSHIFT}_$MissConf.job

MissConfName=`$FindMissConf $Tshift $Ens $jobN | grep "/bgdata4/scnfqcd/miyamoto/conf/pacs-cs"`

count1=0
count2=0
if [ ! -e $gfile ]; then
    while read line; do
	if [ $line != $MissConfName ]; then
            echo $line >> $gfile
	    count1=$((count1+1))
        fi
	count2=$((count2+1))
    done < `cat $BASE1/Tshift${TSHIFT}_${jobN}.job | grep arguments | cut -d " " -f 4`
else
    count2=1
fi

if [ $count1 -ne $((count2-1)) ]; then
    echo "ERROR: Failed to make Gfile. $count1 $count2 $gfile"
    exit 1
fi

if [ ! -e $gfile2 ]; then
    echo $MissConfName > $gfile2
fi

echo "Next,"
echo "(1) bllsubmit $BASE1/Tshift${TSHIFT}_${jobN}_avoid.$MissConf.job"
echo "(2) bllsubmit $BASE1/Tshift${TSHIFT}_$MissConf.job"
