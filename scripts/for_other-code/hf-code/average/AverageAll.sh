#!/bin/bash

if [ $# -ne 5 ]; then
    echo " usage : sh `basename $0` [conf list] [#.average] [Output Tshift (A??, etc.)] [Input directory] [Output directory]"
    exit 1
fi

conf_list=$1
Naverage=$2
oTshift=$3
ibase=$4
obase=$5

######### ! CHANGE HERE ! #########
ScriptsBase=$HOME/analysis_code-set/scripts/for_other-code/hf-code/average

MakeHFdir=$ScriptsBase/../MakeResultsDir-HF.sh

AveCorr=$ScriptsBase/AverageCorr.sh
AveWave=$ScriptsBase/AverageWave.sh

InputDir_List=$ScriptsBase/InputDir.list

CorrPS_SS_List=$ScriptsBase/Corr.PS.SS.list
CorrMulti_List=$ScriptsBase/Corr.multi.list
BBwave_List=$ScriptsBase/BBwave.list
PPwave_List=$ScriptsBase/PPwave.list
DDwave_List=$ScriptsBase/DDwave.list
DOwave_List=$ScriptsBase/DOwave.list

log_dir=$obase/zzz.ave.log
################################### 

mkdir -p $log_dir

$MakeHFdir $obase $conf_list correlator.PS.dir correlator.SS.dir correlator.multi.dir \
    `cat ${BBwave_List}1` `cat ${BBwave_List}2` `cat ${BBwave_List}3` \
    `cat $PPwave_List` `cat $DDwave_List` `cat $DOwave_List`

$AveCorr $conf_list $Naverage $ibase $InputDir_List $obase $oTshift $CorrPS_SS_List 1 > $log_dir/Ave_Corr.PS.SS.log 2>&1 &
$AveCorr $conf_list $Naverage $ibase $InputDir_List $obase $oTshift $CorrMulti_List 2 > $log_dir/Ave_Corr.multi.log 2>&1 &

$AveWave $conf_list $Naverage $ibase $InputDir_List $obase $oTshift ${BBwave_List}1 1 > $log_dir/Ave_BBwave.log1 2>&1 &
$AveWave $conf_list $Naverage $ibase $InputDir_List $obase $oTshift ${BBwave_List}2 1 > $log_dir/Ave_BBwave.log2 2>&1 &
$AveWave $conf_list $Naverage $ibase $InputDir_List $obase $oTshift ${BBwave_List}3 1 > $log_dir/Ave_BBwave.log3 2>&1 &
$AveWave $conf_list $Naverage $ibase $InputDir_List $obase $oTshift $PPwave_List 2    > $log_dir/Ave_PPwave.log  2>&1 &
$AveWave $conf_list $Naverage $ibase $InputDir_List $obase $oTshift $DDwave_List 3    > $log_dir/Ave_DDwave.log  2>&1 &
$AveWave $conf_list $Naverage $ibase $InputDir_List $obase $oTshift $DOwave_List 4    > $log_dir/Ave_DOwave.log  2>&1 &
