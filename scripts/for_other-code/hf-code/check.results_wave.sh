#!/bin/sh

### Auther: Takaya Miyamoto 
### Date  : Fri Sep  2 00:51:57 JST 2016
###
### Brief : To check the file size of the hadron correlator
###       : and show the invalid file size correlator
###

if [ $# -eq 0 ]; then
    echo "usage: `basename $0` [results dir] ..."
    exit 1
fi

ReadBin=$HOME/analysis_code-set/scripts/convenience/Binary_read.py
Length=969

Time="012"

SkipLen=$((4 * (16+$Length)))
for DirName in $@; do
    for BBdir in `ls $DirName | grep BB`; do
	if [ ! -e $DirName/$BBdir ]; then
	    echo "ERROR: $DirName/$BBdir is not exist."
	    exit 1
	fi
	for conf in `ls $DirName/$BBdir`; do
	    for TIME in $Time; do
		NBS=NBSwave.+$TIME+???.000.000.000.$conf.NUC_CG05.NUC_CG05
		tmp=`$ReadBin $DirName/$BBdir/$conf/$NBS $SkipLen double | grep double | awk '{print $3}'`
		if [ `echo $tmp | grep nan | wc -l` -ne 0 ]; then
		    echo "$DirName:$BBdir:$conf:$TIME"
		fi
	    done
	done
    done
done
