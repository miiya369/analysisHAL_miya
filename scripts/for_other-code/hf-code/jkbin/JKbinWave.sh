#!/bin/sh

if [ $# -ne 8 ]; then
    echo " usage : `basename $0` [Org conf list] [JKbin conf list] [Bin size] [NBS list] [TXYZshift (???.???.???.???)] [input dir] [output dir] [1=BBwave, 2=PPwave, 3=DDwave, 4=DOwave]"
    exit 1
fi

######### ! CHANGE HERE ! #########
t_min_BB=6
t_max_BB=18
t_min_PP=6
t_max_PP=31

WaveJKbin_NoComp=$HOME/analysis_code-set/main/indep_code/WaveJKbin_not-compress
WaveJKbin_Comp=$HOME/analysis_code-set/main/indep_code/WaveJKbin_compress
###################################

org_conf_list=$1
JK_conf_list=$2
BinSize=$3
wave_list=$4
TXYZshift=$5
ibase=$6
obase=$7

if [ $8 -eq 1 ]; then
    t_min=$t_min_BB
    t_max=$t_max_BB
    src_rela=CG05
    snk_rela=CG05
    NBStype=NUC
    WaveJKbin=$WaveJKbin_Comp
elif [ $8 -eq 2 ]; then
    t_min=$t_min_PP
    t_max=$t_max_PP
    src_rela=CG05
    snk_rela=CG05
    NBStype=PP
    WaveJKbin=$WaveJKbin_NoComp
elif [ $8 -eq 3 ]; then
    t_min=$t_min_BB
    t_max=$t_max_BB
    src_rela=NR
    snk_rela=NR
    NBStype=DD
    WaveJKbin=$WaveJKbin_NoComp
elif [ $8 -eq 4 ]; then
    t_min=$t_min_BB
    t_max=$t_max_BB
    src_rela=NR
    snk_rela=NR
    NBStype=DO
    WaveJKbin=$WaveJKbin_NoComp
else
    echo "Invalid arguments: Please put 1 - 4 (1=BBwave, 2=PPwave, 3=DDwave, 4=DOwave)"
    exit 1
fi

echo "Setup:"
echo " - Current    = `pwd`"
echo " - Conf List  = $org_conf_list"
echo " - Wave List  = $wave_list"
echo " - TXYZ shift = $TXYZshift"
echo " - Input dir  = $ibase"
echo " - Output dir = $obase"
echo " - time       = $t_min - $t_max"
echo " - Bin size   = $BinSize"
echo
echo "Make JK-bin samples NBSwave START: `date`"
echo

WAVE_LIST=`cat $wave_list`

for WaveDir in $WAVE_LIST
do
    for Time in `seq $t_min 1 $t_max`
    do
	TIME=`printf %03d $Time`
        for SRC_RELA in $src_rela
        do
            for SNK_RELA in $snk_rela
            do		
		ifname=$ibase/$WaveDir/REPCONF/NBSwave.+$TIME+$TXYZshift.REPCONF.${NBStype}_$SNK_RELA.${NBStype}_$SRC_RELA
		ofname=$obase/$WaveDir/REPCONF/NBSwave.+$TIME+$TXYZshift.REPCONF.${NBStype}_$SNK_RELA.${NBStype}_$SRC_RELA
		
		$WaveJKbin $org_conf_list $JK_conf_list $BinSize $ifname $ofname
		if [ $? -ne 0 ]; then
		    echo "ERROR: Something happened. Exit."
		    echo ""
		    echo "Report: time  = $Time"
                    echo "      : dir   = $WaveDir"
                    echo "      : rela  = ${SNK_RELA}_${SRC_RELA}"
		    echo "      : Corr  = $CorrName"
		    echo "      : ifile = $ifname"
		    echo "      : ofile = $ofname"
		    exit 1
		fi
	    done
	done
	
    done
    echo "End: $WaveDir: `date`"
done

echo
echo "Make JK-bin samples NBSwave  END : `date`"
