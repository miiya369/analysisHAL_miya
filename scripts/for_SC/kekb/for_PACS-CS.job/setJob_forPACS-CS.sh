#!/bin/sh

if test $# -ne 3; then
    echo "usage : sh $0 [Tshift] [Ens] [Nconf/job]"
    exit 1
fi

Tshift=$1
Ens=$2
Nconf_job=$3

TSHIFT=`printf %03d $Tshift`

WorkingDir=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF
BridgeDir=$WorkingDir/bridge-1.3.2-hal.v01/src2
makeJobScript=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/script/makeJobscript_forPACS-CS.sh
makeParam=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/script/makeParam.yaml.py
LinkToResults=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/script/LinkToResults.sh

if [ ! -e $WorkingDir/pacs-cs_ens$Ens.RHQ/c_quark/Tshift$TSHIFT ]
then
    mkdir -p $WorkingDir/pacs-cs_ens$Ens.RHQ/c_quark/Tshift$TSHIFT
fi

if [ ! -e $WorkingDir/pacs-cs_ens$Ens.RHQ/s_quark/Tshift$TSHIFT ]
then
    mkdir -p $WorkingDir/pacs-cs_ens$Ens.RHQ/s_quark/Tshift$TSHIFT
fi

if [ ! -e $WorkingDir/pacs-cs_ens$Ens.RHQ/Output/Tshift$TSHIFT ]
then
    mkdir -p $WorkingDir/pacs-cs_ens$Ens.RHQ/Output/Tshift$TSHIFT
fi

if [ ! -e $BridgeDir/pacs-cs_ens$Ens.work/Tshift$TSHIFT ]
then
    cp -r $BridgeDir/pacs-cs_ens$Ens.work/Tshift000 \
	$BridgeDir/pacs-cs_ens$Ens.work/Tshift$TSHIFT
fi

current=`pwd`
cd $BridgeDir/pacs-cs_ens$Ens.work/Tshift$TSHIFT
$LinkToResults $Ens $Tshift
cd $current

$makeJobScript $Tshift $Ens $Nconf_job
$makeParam $Tshift $Ens > $WorkingDir/pacs-cs_ens$Ens.RHQ/Param/quark_Tshift$TSHIFT.yaml

echo "(Next) Submit jobs in pacs-cs_ens$Ens.RHQ/Jobscript/Tshift$TSHIFT"
