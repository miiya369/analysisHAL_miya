#!/bin/sh

if [ $# -ne 4 ]; then
    echo "usage: `basename $0` [Results.tar Dir] [Ens] [0=Clover, 1=RHQ] [Tmod4=?]"
    exit 1
fi

CpUberFTP=$HOME/analysis_code-set/scripts/convenience/cp-uberftp.to.sh

Dir=$1
Ens=$2
fopr=$3
Tmod4=$4

if [ $Ens -gt 3 -o $Ens -lt 1 ]; then
    echo "ERROR: Ens must be put [1-3]"
    exit 1
fi

if [ $fopr -eq 0 ]; then
    if [ $Ens -eq 2 -o $Ens -eq 3 ]; then
	echo "ERROR: You cannot select the Clover when the Ens is not 1."
	exit 1
    fi
fi

if [ $fopr -eq 0 ]; then
    Fopr=Clover
elif [ $fopr -eq 1 ]; then
    Fopr=RHQ
else
    echo "ERROR: Invalid Fopr. Must be put 0=Clover or 1=RHQ."
    exit 1
fi

if [ ! -e $Dir ]; then
    echo "No such a directory '$Dir'"
    exit 1
fi

BaseJLDG=gfarm/pacscs/BBILQCD/miyamoto/results.kekb/results_kekb_pacs-cs_bridge_$Fopr/Ensemble$Ens/TshiftMod04_$Tmod4
#BaseJLDG=gfarm/pacscs/BBILQCD/miyamoto/results.kekb/results_kekb_pacs-cs_bridge_$Fopr/Ensemble$Ens/TshiftMod0x_Average

$CpUberFTP $Dir $BaseJLDG > $Dir.CopyUberftp.to.log 2>&1
