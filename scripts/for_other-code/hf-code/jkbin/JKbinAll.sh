#!/bin/sh

if [ $# -ne 5 ]; then
    echo " usage : sh `basename $0` [Org conf list] [base name of confs] [Bin size] [input dir] [output dir]"
    exit 1
fi

######### ! CHANGE HERE ! #########
TXYZshift=A64.000.000.000

ScriptsBase=$HOME/analysis_code-set/scripts/for_other-code/hf-code/jkbin

MakeJKdir=$ScriptsBase/MakeJKbinDir-HF.sh

JKbinCorr=$ScriptsBase/JKbinCorr.sh
JKbinWave=$ScriptsBase/JKbinWave.sh

Corr_PS_SS_list=$ScriptsBase/Corr.PS.SS.list
Corr_multi_list=$ScriptsBase/Corr.multi.list
BBwave_list=$ScriptsBase/BBwave.list
PPwave_list=$ScriptsBase/PPwave.list
DDwave_list=$ScriptsBase/DDwave.list
DOwave_list=$ScriptsBase/DOwave.list
###################################

org_conf_list=$1
jk_conf_name=$2
BinSize=$3
idir=$4
odir=$5

Nconf_org=`cat $org_conf_list | wc -l`
if [ $? -ne 0 ]; then
    exit 1
fi

echo "Nconf_org = $Nconf_org"
echo "Bin Size  = $BinSize"
echo "#.bin     = $((Nconf_org/BinSize))"

if [ ! -e $odir/zzz.jkbin.log ]; then
    mkdir -p $odir/zzz.jkbin.log
fi 

$MakeJKdir $odir $jk_conf_name $Nconf_org $BinSize `cat $BBwave_list` \
    correlator.PS.dir correlator.SS.dir correlator.multi.dir \
    `cat $PPwave_list` `cat $DDwave_list` `cat $DOwave_list`

jk_conf_list=$ScriptsBase/tmpName_JKbinConf_`LANG=C date +'%H.%M.%S_%b%Y.JST'`
ls $odir/correlator.PS.dir > $jk_conf_list

$JKbinCorr $org_conf_list $jk_conf_list $BinSize $Corr_PS_SS_list $TXYZshift $idir $odir 1 \
    > $odir/zzz.jkbin.log/JKbin_Corr.PS.SS.log 2>&1 &
$JKbinCorr $org_conf_list $jk_conf_list $BinSize $Corr_multi_list $TXYZshift $idir $odir 2 \
    > $odir/zzz.jkbin.log/JKbin_Corr.multi.log 2>&1 &

$JKbinWave $org_conf_list $jk_conf_list $BinSize $BBwave_list $TXYZshift $idir $odir 1 \
    > $odir/zzz.jkbin.log/JKbin_BBwave.log 2>&1 &
$JKbinWave $org_conf_list $jk_conf_list $BinSize $PPwave_list $TXYZshift $idir $odir 2 \
    > $odir/zzz.jkbin.log/JKbin_PPwave.log 2>&1 &
$JKbinWave $org_conf_list $jk_conf_list $BinSize $DDwave_list $TXYZshift $idir $odir 3 \
    > $odir/zzz.jkbin.log/JKbin_DDwave.log 2>&1 &
$JKbinWave $org_conf_list $jk_conf_list $BinSize $DOwave_list $TXYZshift $idir $odir 4 \
    > $odir/zzz.jkbin.log/JKbin_DOwave.log 2>&1 &

# rm $jk_conf_list
