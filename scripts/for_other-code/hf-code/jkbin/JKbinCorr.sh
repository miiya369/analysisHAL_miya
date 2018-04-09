#!/bin/sh

if [ $# -ne 8 ]; then
    echo " usage : `basename $0` [Org conf list] [JKbin conf list] [Bin size] [Corr list] [TXYZshift (???.???.???.???)] [input dir] [output dir] [1=PS/SS, 2=multi]"
    exit 1
fi

######### ! CHANGE HERE ! #########
CorrJKbin=$HOME/analysis_code-set/main/indep_code/CorrJKbin
###################################

org_conf_list=$1
JK_conf_list=$2
BinSize=$3
corr_list=$4
TXYZshift=$5
ibase=$6
obase=$7

if [ $8 -eq 1 ]; then
    CorrDir="PS SS"
elif [ $8 -eq 2 ]; then
    CorrDir=multi
else
    echo "Invalid arguments: Please put 1 or 2 (1=PS/SS, 2=multi)"
    exit 1
fi

echo "Setup:"
echo " - Current    = `pwd`"
echo " - Conf List  = $org_conf_list"
echo " - Corr List  = $corr_list"
echo " - TXYZ shift = $TXYZshift"
echo " - Input dir  = $ibase"
echo " - Output dir = $obase"
echo " - Bin size   = $BinSize"
echo
echo "Make JK-bin samples Correlator START: `date`"
echo

CORR_LIST=`cat $corr_list`

for Dir in $CorrDir
do
    for CorrName in $CORR_LIST
    do	
	ifname=$ibase/correlator.$Dir.dir/REPCONF/${CorrName}.+$TXYZshift.REPCONF
	ofname=$obase/correlator.$Dir.dir/REPCONF/${CorrName}.+$TXYZshift.REPCONF
	
	$CorrJKbin $org_conf_list $JK_conf_list $BinSize $ifname $ofname
	if [ $? -ne 0 ]; then
	    echo "ERROR: Something happened. Exit."
	    echo ""
	    echo "Report: dir   = $Dir"
	    echo "      : Corr  = $CorrName"
	    echo "      : ifile = $ifname"
	    echo "      : ofile = $ofname"
	    exit 1
	fi
    done
    echo "End: correlator.$Dir.dir: `date`"
done

echo
echo "Make JK-bin samples Correlator  END : `date`"
