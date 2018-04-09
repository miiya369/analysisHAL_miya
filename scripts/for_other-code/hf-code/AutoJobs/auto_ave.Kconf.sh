#!/bin/sh

### Auther: Takaya Miyamoto
### Date  : Tue Sep 13 21:09:40 JST 2016
###
### Brief : (1) copy tar file from JLDG
###       : (2) check the md5sum with one on the JLDG
###       : (3) expand tar file and chack files
###       : (4) average
###       : (5) delete the file
###

if [ $# -ne 3 ]; then
    echo "usage: `basename $0` [run-? (h or i)] [traj%10=?? (0 or 5)] [Wave Dir]"
    exit 1
fi

run=$1
traj=$2
WaveDir=$3
Traj=`printf %02d $traj`

# Path should be the ABSOLUTE PATH.
WaveAve=$HOME/analysis_code-set/main/indep_code/WaveAve_compress

base=/work/HALQCD/tmiyamoto/working.dir/tmp/Kconf_results_S1

iDir=$base/run-$run.traj_mod10_$Traj
oDir=$base/average

WaveTslice_min=3
WaveTslice_max=22

TailFname=NUC_CG05.NUC_CG05
Tmod_oFile=Ave

### Arguments check
logDir=$oDir/zzz.logs

case $run in
    h ) 
	case $traj in
	    0 ) Mod=04; Tmods=(0 2); Nave=2; AveIndex=B;;
	    5 ) Mod=04; Tmods=(1 3); Nave=2; AveIndex=B;;
	    * ) echo "Invalid argument: traj%10=$traj"; exit 1;;
	esac ;;
    i ) 
	case $traj in
            0 ) Mod=04; Tmods=(0 2); Nave=2; AveIndex=B;;
            5 ) Mod=24; Tmods=(`seq 1 2 23`); Nave=12; AveIndex=A;;
            * ) echo "Invalid argument: traj%10=$traj"; exit 1;;
        esac ;;
    * ) echo "Invalid argument: run-$run"; exit 1;;
esac

tmp_t=`printf %03d ${Tmods[0]}`
Confs_tar=`ls $iDir/t_mod${Mod}_$tmp_t/$WaveDir`
if [ $? -ne 0 ]; then
    exit 1
fi

if [ ! -d $oDir ]; then
    mkdir -p $oDir
fi
if [ ! -d $logDir ]; then
    mkdir -p $logDir
fi

Nwave=$(($WaveTslice_max - $WaveTslice_min + 1))

exit_code=0

echo
echo "Setup: run-$run, traj%10 = $Traj, WaveDir = $WaveDir"
echo
echo "@@@ JOB start : `date`"
echo

### (1) copy from JLDG
echo -n "@   gfpcopy start "
for Tmod in ${Tmods[@]}; do
    TMOD=`printf %03d $Tmod`
    mkdir -p $oDir/run-$run.traj_mod10_$Traj.t_mod${Mod}_$TMOD
    gfpcopy -p $iDir/t_mod${Mod}_$TMOD/$WaveDir $oDir/run-$run.traj_mod10_$Traj.t_mod${Mod}_$TMOD \
	> $logDir/run-$run.traj_mod10_$Traj.t_mod${Mod}_$TMOD.$WaveDir.gfocopy.log 2>&1
    
    if [ $? -ne 0 ]; then
	echo
        exit 1
    fi
done
echo "- DONE: `date`"

### (2) check md5sum
echo -n "@   md5sum  start "
for Tmod in ${Tmods[@]}; do
    TMOD=`printf %03d $Tmod`
    cd $oDir/run-$run.traj_mod10_$Traj.t_mod${Mod}_$TMOD
    
    cat $iDir/t_mod${Mod}_$TMOD/md5.tar.log | grep $WaveDir > md5.tar.log.jldg
    for conf in $WaveDir/*.tar; do
	md5sum ./$conf >> md5.tar.log.here
    done
    
    if [ `cat md5.tar.log.jldg | wc -l` -ne `cat md5.tar.log.here | wc -l` ]; then
	echo
	echo "Invalid #.conf: t_mod${MOD}_$TMOD"
	exit 1
    fi
    if [ `diff md5.tar.log.???? | wc -w` -ne 0 ]; then
	echo
	echo "Invalid check sum: t_mod${MOD}_$TMOD"
	exit 1
    fi
    rm md5.tar.log.????
done
echo "- DONE: `date`"

### (3) expand tar file
echo -n "@   tar     start "
for Tmod in ${Tmods[@]}; do
    TMOD=`printf %03d $Tmod`
    cd $oDir/run-$run.traj_mod10_$Traj.t_mod${Mod}_$TMOD/$WaveDir
    
    for conf in *; do
	tar xf $conf
	fname=`basename $conf .tar`
	if [ `ls $fname | wc -l` -ne $Nwave ]; then
            echo
	    echo "Invalid #.wave (`ls $fname | wc -l` != $Nwave): t_mod${MOD}_$TMOD, $fname"
            exit 1
	fi
	rm $conf
    done
done
echo "- DONE: `date`"

### (4) average
echo -n "@   average start "
if [ ! -d $oDir/results.ave/$WaveDir ]; then
    mkdir -p $oDir/results.ave/$WaveDir
fi
cd $oDir
 
for conf_tar in $Confs_tar; do
    conf=`basename $conf_tar .tar`
    if [ ! -d results.ave/$WaveDir/$conf ]; then
	mkdir -p results.ave/$WaveDir/$conf
    fi
    for time in `seq $WaveTslice_min 1 $WaveTslice_max`; do
	TIME=`printf %03d $time`
	
	iFname=""
	iFnameBase1=run-$run.traj_mod10_$Traj.t_mod${Mod}_
	iFnameBase2=/$WaveDir/$conf/NBSwave.+$TIME+$AveIndex
	iFnameBase3=.000.000.000.$conf.$TailFname
	for Tmod in ${Tmods[@]}; do
	    TMOD2=`printf %02d $Tmod`
	    TMOD3=`printf %03d $Tmod`
	    iFname="$iFname $iFnameBase1$TMOD3$iFnameBase2$TMOD2$iFnameBase3"
	done
	
	if [ `echo $iFname | wc -w` -ne $Nave ]; then
	    echo
	    echo "Invalid #.ave (`echo $iFname | wc -w` != $Nave): $conf, t=$time"
	    exit 1
	fi
	
	oFname=results.ave/$WaveDir/$conf/NBSwave.+$TIME+$Tmod_oFile.000.000.000.$conf.$TailFname
	
	$WaveAve $oFname $iFname
	if [ $? -ne 0 ]; then
	    echo
            echo "Error while averaging."
	    exit 1
	fi
    done
done
echo "- DONE: `date`"

### (5) delete files
echo -n "@   delete  start "
for Tmod in ${Tmods[@]}; do
    TMOD=`printf %03d $Tmod`
    rm -r run-$run.traj_mod10_$Traj.t_mod${Mod}_$TMOD/$WaveDir
done
echo "- DONE: `date`"

echo
echo "@@@ JOB  end  : `date`"
echo

exit $exit_code
