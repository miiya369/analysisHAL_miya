#!/bin/sh

if [ $# -ne 3 ]
then
    echo "usage : sh $0 [#.conf] [#.conf/job] [Job Num.]"
    exit 1
fi

Nconf=$1
Nconf_job=$2
JobN=$3

######### ! CHANGE HERE ! #########
JobName=HF-CP_$JobN
WorkDir=/bghome/scnfqcd/tkymiya/tmp/for_cp-pacs/bridge-1.3.2-hal.v01/src2/run.cp-pacs
Tshift=0
que=32
jobc=("a" "b" "c")
###################################

TSHIFT=`printf %03d $Tshift`

Njobc=${#jobc[@]}

conf_s=$((JobN   * Nconf_job + 1))
conf_e=$((conf_s + Nconf_job - 1))
if [ $conf_e -gt $Nconf ]
then
    conf_e=$Nconf
fi

conf=`printf %03d $conf_s`-`printf %03d $conf_e`

echo "jobclass=qb$que${jobc[$((JobN % Njobc))]}"
echo "inputfile="
echo "workingdir=$WorkDir"
echo "jobname=$JobName"
echo "outputfile=\${workingdir}/Output/run.out_$JobN"
echo "errorfile=\${workingdir}/Output/run.err_$JobN"
echo "executable=$WorkDir/run.x"
echo 'environment="OMP_NUM_THREADS=8:MUSPI_NUMINJFIFOS=8:MUSPI_NUMBATIDS=10"'
echo "tasks_per_node=8"
echo "mapfile="
echo "arguments=-params $WorkDir/input.yaml/quark.yaml -f $WorkDir/gfile.lst/gfile_list.pacs-cs+JLQCD_ens3.gfix_$conf"
