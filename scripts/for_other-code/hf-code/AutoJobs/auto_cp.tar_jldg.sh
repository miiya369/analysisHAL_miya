#!/bin/sh

### Auther: Takaya Miyamoto 
### Date  : Tue Sep 13 21:09:40 JST 2016
###
### Brief : (1) copy tar file from JLDG
###       : (2) check the md5sum with one on the JLDG
###       : (3) expand tar file and check files
###       : (4) tar c (for HAL) and check tar files
###       : (5) copy new tar files to JLDG (Replace old tar file)
###       : (6) delete the file (Back Ground)


if [ $# -ne 4 ]; then
    echo "usage: `basename $0` [Ens] [s/c] [Tmod] [Tsrc (If run all time, set '999')]"
    exit 1
fi

TarForHAL=$HOME/analysis_code-set/scripts/for_other-code/hf-code/tar_c.hal.sh
CheckResults=$HOME/analysis_code-set/scripts/for_other-code/hf-code/check.results.sh
CheckTarHAL=$HOME/analysis_code-set/scripts/for_other-code/hf-code/check.tar.hal.sh
CheckTarLogHAL=$HOME/analysis_code-set/scripts/for_other-code/hf-code/check.tar.log.hal.sh

Ens=$1
Qrk=$2
Tmod=$3
tsrc=$4

### Arguments check
if [ $Ens -lt 1 -o $Ens -gt 3 ]; then
    echo "Arguments error: You should set Ens = [1-3]"
    exit 1
elif [ $Qrk != "s" -a $Qrk != "c" ]; then
    echo "Arguments error: You should set [s, c] for quark"
    exit 1
elif [ $Tmod -lt 0 -o $Tmod -gt 3 ]; then
    echo "Arguments error: You should set Tmod = [0-3]"
    exit 1
elif [ $Ens -ne 1 ]; then
    if [ $Qrk = "s" ]; then
	echo "Arguments error: You don't set Ens != 1 and quark = s"
	exit 1
    fi
fi

case $Ens in
    1) LinesTarLog=1507884; SizeTarLog=271016543; NumFiles=1483482;;
    2) LinesTarLog=1511663; SizeTarLog=271695778; NumFiles=1487200;;
    3) LinesTarLog=1700613; SizeTarLog=305657528; NumFiles=1673100;;
esac
case $Qrk in
    "s") Fopr=Clover;;
    "c") Fopr=RHQ;;
esac
Base=/tmp/miyamoto.jldg/disk.jldg/gfarm/pacscs/BBILQCD/miyamoto/results.kekb/results_kekb_pacs-cs_bridge_$Fopr/Ensemble$Ens/TshiftMod04_$Tmod

CheckTarHAL=$CheckTarHAL.ens$Ens
CheckTarLogHAL=$CheckTarLogHAL.ens$Ens

if [ $tsrc -eq 999 ]; then
    tsrc=`seq $Tmod 4 63`
fi

exit_code=0

echo "@@@ JOB start : `date`"
echo

for Tsrc in $tsrc; do
    TSRC=`printf %03d $Tsrc`
    echo "@@  Tsrc = $TSRC start : `date`"

### (1) copy from JLDG
    echo -n "@   cp from start : `date`"
    if [ ! -f Tshift$TSRC.tar -a ! -d Tshift$TSRC ]; then
	gfpcopy -p $Base/Tshift$TSRC.tar . > Tshift$TSRC.tar.gfpcopy.from.log 2>&1
	if [ $? -ne 0 ]; then
	    echo
	    echo "ERROR while gfpcopy from"
	    echo
	    exit_code=1
            continue
	fi
	rm Tshift$TSRC.tar.gfpcopy.from.log
    fi
    echo " - Successful"

### (2) check md5sum
    echo -n "@   md5sum  start : `date`"
    if [ -f Tshift$TSRC.tar ]; then
	cksum_here=`md5sum Tshift$TSRC.tar | awk '{print $1}'`
	cksum_jldg=`gfcksum $Base/Tshift$TSRC.tar | awk '{print $1}'`
	if [ $cksum_here != $cksum_jldg ]; then
	    echo
	    echo "ERROR while md5sum: Differect md5sum"
	    echo
	    exit_code=1
            continue
	fi
    fi
    echo " - Successful"

### (3)-1 expand tar file
    echo -n "@   tar x   start : `date`"
    if [ -f Tshift$TSRC.tar -a ! -d Tshift$TSRC ]; then
	tar xvf Tshift$TSRC.tar > Tshift$TSRC.tar.tar.log 2>&1
	if [ $? -ne 0 ]; then
	    echo
            echo "ERROR while tar x: tar x is not successful"
            echo
	    exit_code=1
            continue
	fi
	if [ `wc -l Tshift$TSRC.tar.tar.log | awk '{print $1}'` -ne $LinesTarLog ]; then
	    echo
	    echo "ERROR while tar x: different #.lines of tar.log != $LinesTarLog"
	    echo
	    exit_code=1
            continue
	fi
	if [ `wc -c Tshift$TSRC.tar.tar.log | awk '{print $1}'` -ne $SizeTarLog ]; then
	    echo
	    echo "ERROR while tar x: different size of tar.log != $SizeTarLog"
	    echo
	    exit_code=1
            continue
	fi
	rm Tshift$TSRC.tar
	rm Tshift$TSRC.tar.tar.log
    fi
    echo " - Successful"

### (3)-2 Check results file
    echo -n "@   chk res start : `date`"
    if [ -d Tshift$TSRC ]; then
	if [ `$CheckResults Tshift$TSRC | wc -l` -ne 0 ]; then
	    echo
	    echo "ERROR while checking the results/: There are invalid file(s)"
	    echo
	    exit_code=1
            continue
	fi
	if [ `find Tshift$TSRC -type f | wc -l` -ne $NumFiles ]; then
            echo
            echo "ERROR while checking the results/: Invalid #.file"
            echo
	    exit_code=1
            continue
	fi
    fi
    echo " - Successful"
    
### (4)-1 compress tar file (for HAL)
    echo -n "@   tar c   start : `date`"
    if [ ! -d Tshift$TSRC.tar ]; then
	$TarForHAL Tshift$TSRC
	if [ $? -ne 0 ]; then
	    echo
            echo "ERROR while tar c"
            echo
	    exit_code=1
            continue
	fi
    fi
    echo " - Successful"

### (4)-2 Check results tar file (for HAL)
    echo -n "@   chk tar start : `date`"
    if [ -d Tshift$TSRC.tar ]; then
	if [ `ls Tshift$TSRC.tar/*.tar | wc -l` -ne 59 -o `ls Tshift$TSRC.tar/zzz.tar.log | wc -l` -ne 59 ]; then
            echo
            echo "ERROR while checking the results.tar/: Invalid #.tar files"
            echo
            exit_code=1
            continue
	fi
	if [ `$CheckTarHAL Tshift$TSRC.tar | grep -v correct | wc -l` -ne 0 ]; then
            echo
            echo "ERROR while checking the results.tar/: The tar file(s) have invalid size"
            echo
	    exit_code=1
            continue
	fi
	if [ `$CheckTarLogHAL Tshift$TSRC.tar/zzz.tar.log | grep -v correct | wc -l` -ne 0 ]; then
            echo
            echo "ERROR while checking the results.tar/zzz.tar.log: The tar.log have invalid #.line"
            echo
	    exit_code=1
            continue
	fi
    fi
    echo " - Successful"

### (5) copy to jldg
    echo -n "@   cp to   start : `date`"
    if [ -d Tshift$TSRC.tar ]; then
	gfmv $Base/Tshift$TSRC.tar $Base/Tshift$TSRC.tar.org
	gfpcopy -p Tshift$TSRC.tar $Base > Tshift$TSRC.tar.gfpcopy.to.log 2>&1
	if [ $? -ne 0 ]; then
	    echo
	    echo "ERROR while gfpcopy to"
	    echo
	    exit_code=1
	    continue
	fi
	rm Tshift$TSRC.tar.gfpcopy.to.log
    fi
    echo " - Successful"

### (6) delete files
    echo "@   delete  start : Back Ground"
    rm -r Tshift$TSRC.tar &
    rm -r Tshift$TSRC &

    echo "@@  Tsrc = $TSRC  end  : `date`"
    echo
done

echo "@@@ JOB  end  : `date`"

exit $exit_code
