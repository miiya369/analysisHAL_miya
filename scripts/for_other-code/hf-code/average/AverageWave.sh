#!/bin/sh

if [ $# -ne 8 ]; then
    echo " usage : sh `basename $0` [conf list] [#.average] [Input directory base] [Input directory list] [Output directory] [Output Tshift (A??, etc.)] [NBSwave list] [1=BBwave, 2=PPwave, 3=DDwave, 4=DOwave]"
    exit 1
fi

######### ! CHANGE HERE ! #########
t_min_BB=6
t_max_BB=18
t_min_PP=6
t_max_PP=31

XYZshift=000.000.000

WaveAve_NoComp=$HOME/analysis_code-set/main/indep_code/WaveAve_not-compress
WaveAve_Comp=$HOME/analysis_code-set/main/indep_code/WaveAve_compress
################################### 

conf_list=$1
Naverage=$2
ibase=$3
idir_list=$4
obase=$5
oTshift=$6
wave_list=$7

if [ $8 -eq 1 ]; then
    t_min=$t_min_BB
    t_max=$t_max_BB
    src_rela=CG05
    snk_rela=CG05
    NBStype=NUC
    WaveAve=$WaveAve_Comp
elif [ $8 -eq 2 ]; then
    t_min=$t_min_PP
    t_max=$t_max_PP
    src_rela=CG05
    snk_rela=CG05
    NBStype=PP
    WaveAve=$WaveAve_NoComp
elif [ $8 -eq 3 ]; then
    t_min=$t_min_BB
    t_max=$t_max_BB
    src_rela=NR
    snk_rela=NR
    NBStype=DD
    WaveAve=$WaveAve_NoComp
elif [ $8 -eq 4 ]; then
    t_min=$t_min_BB
    t_max=$t_max_BB
    src_rela=NR
    snk_rela=NR
    NBStype=DO
    WaveAve=$WaveAve_NoComp
else
    echo "Invalid arguments: Please put 1 - 4 (1=BBwave, 2=PPwave, 3=DDwave, 4=DOwave)"
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
if [ ! -f $wave_list ]; then
    echo "ERROR: $wave_list is not found"
    exit 1
fi

if [ $Naverage -ne `cat $idir_list | wc -l` ]; then
    echo "ERROR: #.Input directory != #.average"
    exit 1
fi

for idir in `cat $idir_list`; do
    for WaveDir in `cat $wave_list`; do
        if [ ! -d $ibase/$idir/$WaveDir ]; then
            echo "ERROR: $ibase/$idir/$WaveDir is not found"
            exit 1
        fi
    done
done

###### Average Start ######
echo "Setup:"
echo " - Current    = `pwd`"
echo " - Conf  List = $conf_list"
echo " - Wave  List = $wave_list"
echo " - Input Base = $ibase"
echo " - Input List = $idir_list"
echo " - Output dir = $obase"
echo " - time       = $t_min - $t_max"
echo " - Out Tshift = $oTshift"
echo " - #.average  = $Naverage"
echo
echo "Average NBSwave START: `date`"
echo

CONF_LIST=`cat $conf_list`
WAVE_LIST=`cat $wave_list`
INPUT_LIST=`cat $idir_list`

for CONF in $CONF_LIST
do
    for Time in `seq $t_min 1 $t_max`
    do
	TIME=`printf %03d $Time`
	for WaveDir in $WAVE_LIST
	do
	    for SRC_RELA in $src_rela
	    do
		for SNK_RELA in $snk_rela
		do
		    ofname=$obase/$WaveDir/$CONF/NBSwave.+$TIME+$oTshift.$XYZshift.$CONF.${NBStype}_$SNK_RELA.${NBStype}_$SRC_RELA
		    ifname=""
		    
		    for idir in $INPUT_LIST; do
			ifname="$ifname `echo $ibase/$idir/$WaveDir/$CONF/NBSwave.+$TIME+???.$XYZshift.$CONF.${NBStype}_$SNK_RELA.${NBStype}_$SRC_RELA`"
		    done
		    
		    if [ `echo $ifname | wc -w` -ne $Naverage ]; then
			echo "ERROR: Invalid #.average, `echo $ifname | wc -w` != $Naverage. Exit."
			echo ""
			echo "Report: conf = $CONF"
                        echo "      : time = $Time"
                        echo "      : dir  = $WaveDir"
                        echo "      : rela = ${SNK_RELA}_${SRC_RELA}"
                        exit 1
                    fi
		    
		    $WaveAve $ofname $ifname
		    if [ $? -ne 0 ]
		    then
			echo "ERROR: Something happened. Exit."
			echo
			echo "Report: conf = $CONF"
			echo "      : time = $Time"
			echo "      : dir  = $WaveDir"
			echo "      : rela = ${SNK_RELA}_${SRC_RELA}"
			exit 1
		    fi
		done
	    done
	done
    done
    echo "End: $CONF: `date`"
done

echo
echo "Average NBSwave  END : `date`"
