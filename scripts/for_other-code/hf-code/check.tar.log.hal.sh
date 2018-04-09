#!/bin/sh

if [ $# -eq 0 ]; then
    echo "usage: `basename $0` [tar.log dir] ..."
    exit 1
fi

######### Setup #########
Nconf=399

#calcT_BB_min=0
#calcT_BB_max=31
#calcT_PP_min=0
#calcT_PP_max=31
#calcT_DD_min=0
#calcT_DD_max=31
#calcT_DO_min=0
#calcT_DO_max=31

calcT_BB_min=6
calcT_BB_max=18
calcT_PP_min=6
calcT_PP_max=31
calcT_DD_min=6
calcT_DD_max=18
calcT_DO_min=6
calcT_DO_max=18

# N.type means how many CG05 and NR in a NBS
Ntype_BB=2
#Ntype_BB=1
Ntype_PP=1
Ntype_DD=1
Ntype_DO=1

Nspin_DD=1
Nspin_DO=2

Ncorr_PS_SS=40
Ncorr_multi=208
#Ncorr_multi=210

Nfile_NoSplit_multi=5
#########################

BBsize=$(($Nconf*$Ntype_BB*($calcT_BB_max-$calcT_BB_min+1)+$Nconf+1))
PPsize=$(($Nconf*$Ntype_PP*($calcT_PP_max-$calcT_PP_min+1)+$Nconf+1))
DDsize=$(($Nconf*$Ntype_DD*($calcT_DD_max-$calcT_DD_min+1)*$Nspin_DD+$Nconf*$Nspin_DD+$Nspin_DD+1))
DOsize=$(($Nconf*$Ntype_DO*($calcT_DO_max-$calcT_DO_min+1)*$Nspin_DO+$Nconf*$Nspin_DO+$Nspin_DO+1))

PS_SSsize=$(($Nconf*$Ncorr_PS_SS+$Nconf+1))
multisize=$(($Nconf*$Ncorr_multi+$Nconf+1))
nosplsize=$(($Nconf*(1+1+$Nfile_NoSplit_multi)+$Nconf*3+3+1))

for DirName in $@; do
    Files=`ls $DirName`
    
    if [ $? -ne 0 ]; then
	exit 1
    fi
    
    count=0
    
    for File in $Files; do
	Fsize=`wc -l $DirName/$File | cut -d " " -f 1`
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
	    if [ $Fsize -ne $PS_SSsize ]; then
		echo "$DirName/$File: $Fsize != $PS_SSsize"
		count=$((count+1))
            fi
	elif [ `echo $File | grep "correlator.SS"` ]; then
            if [ $Fsize -ne $PS_SSsize ]; then
		echo "$DirName/$File: $Fsize != $PS_SSsize"
		count=$((count+1))
            fi
	elif [ `echo $File | grep "correlator.multi"` ]; then
            if [ $Fsize -ne $multisize ]; then
		echo "$DirName/$File: $Fsize != $multisize"
		count=$((count+1))
            fi
	elif [ `echo $File | grep "no_split"` ]; then
            if [ $Fsize -ne $nosplsize ]; then
		echo "$DirName/$File: $Fsize != $nosplsize"
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
