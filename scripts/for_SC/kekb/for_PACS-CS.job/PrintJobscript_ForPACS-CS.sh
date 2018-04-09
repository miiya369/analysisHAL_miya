#!/bin/sh

if test $# -ne 6; then
    echo "usage : sh $0 [Tshift] [Job Num] [Job class] [Nconf Start] [Nconf End] [Ens]"
    exit 1
fi

Ens=$6

if [ $Ens == "1" ]
then
    que=32
#    que=128
elif [ $Ens == "2" ]
then
    que=128
elif [ $Ens == "3" ]
then
    que=128
else
    exit 1
fi

Tshift=$1
JobN=$2
TSHIFT=`printf %03d $Tshift`
jobc=$3
conf_s=$4
conf_e=$5

conf=`printf %03d $conf_s`-`printf %03d $conf_e`

echo "jobclass=qb$que$jobc"
echo "inputfile="
echo "workingdir=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/bridge-1.3.2-hal.v01/src2/pacs-cs_ens$Ens.work/Tshift$TSHIFT"
echo "jobname=PC${Ens}_T${Tshift}_$JobN"
echo "outputfile=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/pacs-cs_ens$Ens.RHQ/Output/Tshift$TSHIFT/run.out.Tshift${TSHIFT}_$JobN"
echo "errorfile=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/pacs-cs_ens$Ens.RHQ/Output/Tshift$TSHIFT/run.err.Tshift${TSHIFT}_$JobN"
echo "executable=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/bridge-1.3.2-hal.v01/src2/pacs-cs_ens$Ens.work/Tshift$TSHIFT/run.x"
echo 'environment="OMP_NUM_THREADS=8:MUSPI_NUMINJFIFOS=8:MUSPI_NUMBATIDS=10"'
echo "tasks_per_node=8"
echo "mapfile="
echo "arguments=-params /bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/pacs-cs_ens$Ens.RHQ/Param/quark_Tshift$TSHIFT.yaml -f /bghome/scnfqcd/tkymiya/gfile_list/pacs-cs_ens${Ens}_gfix/gfile_list.pacs-cs_ens${Ens}_gfix_$conf"
