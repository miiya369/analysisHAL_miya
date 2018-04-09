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

size_min=15000
size_max=15600

had_dir="PS"

hadron="pion phi"
corr=Kaon_correlator.+???.000.000.000

for DirName in $@; do
    for HadDir in $had_dir; do
	if [ ! -e $DirName/correlator.$HadDir.dir ]; then
	    echo "ERROR: $DirName/correlator.$HadDir.dir is not exist."
	    exit 1
	fi
	for conf in `ls $DirName/correlator.$HadDir.dir`; do
	    for Hadron in $hadron; do
		corr=${Hadron}_correlator.+???.000.000.000.$conf
		size=`wc -c $DirName/correlator.$HadDir.dir/$conf/$corr | cut -d " " -f 1`
		tmp=`expr $size + 1` > /dev/null 2>&1
		if [ $? -ne 0 -o $size -lt $size_min -o $size -gt $size_max ]; then
		    echo "$DirName:$HadDir:$conf:$Hadron: $size"
		fi
	    done
	done
    done
done
