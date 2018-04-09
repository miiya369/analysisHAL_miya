#!/bin/sh

### Auther: Takaya Miyamoto 
### Date  : Tue Sep 13 21:09:40 JST 2016
###
### Brief : (1) copy tar file from JLDG
###       : (2) check the md5sum with one on the JLDG
###       : (3) expand tar file and chack files
###       : (4) average
###       : (5) delete the file (Back Ground)
###

if [ $# -ne 2 ]; then
    echo "usage: `basename $0` [Tmod (If run all mod, set '999')] [BB.dir (?.??)]"
    exit 1
fi

LinesTarLog=410001
SizeTarLog_not3=61877220
SizeTarLog_3=61057218
NumFiles=409600
NumFiles_ave=10400

MakeHFdir=$HOME/Datas/pacs-cs/ens2/MakeResultsDir-HF.sh
AveWave=$HOME/Datas/pacs-cs/ens2/AverageWave.ens2.mod4.sh
conf_list=$HOME/Datas/pacs-cs/ens2/conf_list.pacs-cs.ens2

tmod=$1
BBdir=$2

### Arguments check
if [ $tmod -lt 0 -o $tmod -gt 3 ]; then
    if [ $tmod -eq 999 ]; then
	tmod="0 1 2 3"
    else
	echo "Arguments error: You should set Tmod = [0-3], or 999 for all mod."
	exit 1
    fi
fi

exit_code=0

echo "@@@ JOB start : `date`"
echo

for Tmod in $tmod; do
    Base1=/tmp/miyamoto.jldg/disk.jldg/gfarm/pacscs/BBILQCD/kesasaki/results/HF.Esb2/mod0$Tmod
    FileN=HF.Esb2.mod0$Tmod.BBwave.dir.S$BBdir.tar
    Base2=/home/LATTICE/student/miyamoto/Datas/pacs-cs/ens2/by_Sasaki-san/TshiftAve/mod4_$Tmod
    
    if [ $Tmod -eq 3 ]; then 
	SizeTarLog=$SizeTarLog_3
    else 
	SizeTarLog=$SizeTarLog_not3
    fi
    
    echo "@@  Tmod$Tmod start : `date`"
    
### (1) copy from JLDG
    echo -n "@   cp from start : `date`"
    if [ ! -f $FileN -a ! -d Tmod$Tmod.S$BBdir ]; then
	gfpcopy -p $Base1/$FileN . > Tmod$Tmod.S$BBdir.gfpcopy.from.log 2>&1
	if [ $? -ne 0 ]; then
	    echo
	    echo "ERROR while gfpcopy from"
	    echo
	    exit_code=1
	    continue
	fi
	rm Tmod$Tmod.S$BBdir.gfpcopy.from.log
    fi
    echo " - Successful"
    
### (2) check md5sum
    echo -n "@   md5sum  start : `date`"
    if [ -f $FileN ]; then
	cksum_here=`md5sum $FileN | awk '{print $1}'`
	cksum_jldg=`gfcksum $Base1/$FileN | awk '{print $1}'`
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
    if [ -f $FileN -a ! -d Tmod$Tmod.S$BBdir ]; then
	mkdir -p Tmod$Tmod.S$BBdir; cd Tmod$Tmod.S$BBdir
	tar xvf ../$FileN > Tmod$Tmod.S$BBdir.tar.log 2>&1
	if [ $? -ne 0 ]; then
	    echo
	    echo "ERROR while tar x: tar x is not successful"
	    echo
	    exit_code=1
	    continue
	fi
	if [ `wc -l Tmod$Tmod.S$BBdir.tar.log | awk '{print $1}'` -ne $LinesTarLog ]; then
	    echo
	    echo "ERROR while tar x: different #.lines of tar.log != $LinesTarLog"
	    echo
	    exit_code=1
	    continue
	fi
	if [ `wc -c Tmod$Tmod.S$BBdir.tar.log | awk '{print $1}'` -ne $SizeTarLog ]; then
	    echo
	    echo "ERROR while tar x: different size of tar.log != $SizeTarLog"
	    echo
	    exit_code=1
	    continue
	fi
	cd ../
	rm $FileN
    fi
    echo " - Successful"

### (3)-2 Check results file
    echo -n "@   chk res start : `date`"
    if [ -d Tmod$Tmod.S$BBdir ]; then
	if [ `find Tmod$Tmod.S$BBdir/BBwave.dir.S$BBdir -type f | wc -l` -ne $NumFiles ]; then
            echo
            echo "ERROR while checking the results/: Invalid #.file"
            echo
	    exit_code=1
            continue
	fi
    fi
    echo " - Successful"
    
### (4) average
    echo -n "@   average start : `date`"
    if [ ! -d $$Base2/org/BBwave.dir.S$BBdir ]; then
	$MakeHFdir $Base2/org $conf_list BBwave.dir.S$BBdir
    fi
    $AveWave Tmod$Tmod.S$BBdir $Base2/org BBwave.dir.S$BBdir $Tmod > $Base2/log/Ave_BBwave.S$BBdir.log 2>&1
    if [ $? -ne 0 ]; then
	echo
	echo "ERROR while average"
	echo
	exit_code=1
	continue
    fi
    if [ `find $Base2/org/BBwave.dir.S$BBdir -type f | wc -l` -ne $NumFiles_ave ]; then
        echo
        echo "ERROR while checking the results.ave/: Invalid #.file"
        echo
	exit_code=1
        continue
    fi
    echo " - Successful"
    
### (5) delete files
    echo "@   delete  start : Back Ground"
    rm -r Tmod$Tmod.S$BBdir &
    
    echo "@@  Tmod$Tmod  end  : `date`"
    echo
done

echo "@@@ JOB  end  : `date`"

exit $exit_code
