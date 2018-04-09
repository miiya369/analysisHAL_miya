#!/bin/sh

if [ $# -eq 0 ]; then
    echo "usage: `basename $0` [tar dir] ..."
    exit 1
fi

######### Setup #########
BBsize=3242670080
PPsize=6713927680
DDsize=6713927680
DOsize=13427855360

corPSsize=260M
corSSsize=262M
multisize=755M
nosplsize=1.1G
#########################

for DirName in $@; do
    Files=`ls $DirName/*.tar`
    
    if [ $? -ne 0 ]; then
	exit 1
    fi
    
    count=0
    
    for File in $Files; do
	Fsize=`ls -l $File | cut -d " " -f 5`
	FSIZE=`ls -lh $File | cut -d " " -f 5`
	if [ `echo $File | grep "BBwave.dir"` ]; then
	    if [ $Fsize -ne $BBsize ]; then
		echo "$DirName/$File: $Fsize != $BBsize"
		count=$((count+1))
	    fi
	elif [ `echo $File | grep "PPwave.dir"` ]; then
            if [ $Fsize -ne $PPsize ]; then
		echo "$DirName/$File: $Fsize != $PPsize"
		count=$((count+1))
            fi
	elif [ `echo $File | grep "Proj.DDwave.dir"` ]; then
            if [ $Fsize -ne $DDsize ]; then
		echo "$DirName/$File: $Fsize != $DDsize"
		count=$((count+1))
            fi
	elif [ `echo $File | grep "Proj.DOwave.dir"` ]; then
            if [ $Fsize -ne $DOsize ]; then
		echo "$DirName/$File: $Fsize != $DOsize"
		count=$((count+1))
            fi
	elif [ `echo $File | grep "correlator.PS"` ]; then
	    if [ $FSIZE != $corPSsize ]; then
		echo "$DirName/$File: $FSIZE != $corPSsize"
		count=$((count+1))
            fi
	elif [ `echo $File | grep "correlator.SS"` ]; then
            if [ $FSIZE != $corSSsize ]; then
		echo "$DirName/$File: $FSIZE != $corSSsize"
		count=$((count+1))
            fi
	elif [ `echo $File | grep "correlator.multi"` ]; then
            if [ $FSIZE != $multisize ]; then
		echo "$DirName/$File: $FSIZE != $multisize"
		count=$((count+1))
            fi
	elif [ `echo $File | grep "no_split"` ]; then
            if [ $FSIZE != $nosplsize ]; then
		echo "$DirName/$File: $FSIZE != $nosplsize"
		count=$((count+1))
            fi
	else
	    echo "Invalid file is included: $DirName/$File"
	    count=$((count+1))
	fi
    done
    
    if [ $count -eq 0 ]; then
	echo "$DirName: correct"
    fi
done
