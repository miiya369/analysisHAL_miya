#!/bin/sh

if [ $# -ne 8 ]; then
    echo " usage : sh `basename $0` [conf list] [#.average] [Input directory base] [Input directory list] [Output directory] [Output Tshift (A??, etc.)] [Corr list] [1=PS/SS, 2=multi]"
    exit 1
fi

######### ! CHANGE HERE ! #########
XYZshift=000.000.000

CorrAve=$HOME/analysis_code-set/scripts/for_other-code/hf-code/average/Scripts/corrave.py
################################### 

conf_list=$1
Naverage=$2
ibase=$3
idir_list=$4
obase=$5
oTshift=$6
corr_list=$7

if [ $8 -eq 1 ]; then
    CorrDir="PS SS"
elif [ $8 -eq 2 ]; then
    CorrDir=multi
else
    echo "Invalid arguments: Please put 1 or 2 (1=PS/SS, 2=multi)"
    exit 1
fi

if [ ! -f $conf_list ]; then ### Error check
    echo "ERROR: $conf_list is not found"
    exit 1
fi
if [ ! -f $idir_list ]; then
    echo "ERROR: $idir_list is not found"
    exit 1
fi
if [ ! -f $corr_list ]; then
    echo "ERROR: $corr_list is not found"
    exit 1
fi

if [ $Naverage -ne `cat $idir_list | wc -l` ]; then
    echo "ERROR: #.Input directory != #.average"
    exit 1
fi

for idir in `cat $idir_list`; do
    for Dir in $CorrDir; do
	if [ ! -d $ibase/$idir/correlator.$Dir.dir ]; then
	    echo "ERROR: $ibase/$idir/correlator.$Dir.dir is not found"
	    exit 1
	fi
    done
done

###### Average Start ######
echo "Setup:"
echo " - Current    = `pwd`"
echo " - Conf  List = $conf_list"
echo " - Corr  List = $corr_list"
echo " - Input Base = $ibase"
echo " - Input List = $idir_list"
echo " - Output dir = $obase"
echo " - Out Tshift = $oTshift"
echo " - #.average  = $Naverage"
echo
echo "Average Correlator START: `date`"
echo

CONF_LIST=`cat $conf_list`
CORR_LIST=`cat $corr_list`
INPUT_LIST=`cat $idir_list`

for CONF in $CONF_LIST
do
    for Dir in $CorrDir
    do
	for CorrName in $CORR_LIST
	do    
	    ofname=$obase/correlator.$Dir.dir/$CONF/${CorrName}.+$oTshift.$XYZshift.$CONF
	    ifname=""
	    
	    for idir in $INPUT_LIST; do
                ifname="$ifname `echo $ibase/$idir/correlator.$Dir.dir/$CONF/${CorrName}.+???.$XYZshift.$CONF`"
            done
	    
	    if [ `echo $ifname | wc -w` -ne $Naverage ]; then
                echo "ERROR: Invalid #.average, `echo $ifname | wc -w` != $Naverage. Exit."
                echo ""
                echo "Report: conf = $CONF"
                echo "      : dir  = $Dir"
                echo "      : Corr = $CorrName"
		exit 1
            fi
	    
	    $CorrAve $ofname $ifname
	    if [ $? -ne 0 ]; then
		echo "ERROR: Something happened. Exit."
		echo ""
		echo "Report: conf = $CONF"
		echo "      : dir  = $Dir"
		echo "      : Corr = $CorrName"
		exit 1
	    fi
	done
    done
    echo "End: $CONF: `date`"
done

echo
echo "Average Correlator  END : `date`"
