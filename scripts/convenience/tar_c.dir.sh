#!/bin/sh

# Author: Takaya Miyamoto
# Date  : Tue Oct 11 20:53:09 JST 2016
#
# Brief : (1) Take tar-compress for each file in the certain directory
#       : (2) Store the log of tar into 'zzz.tar.log'
#       : (3) Store the md5sum value and file size into 'check'

if [ $# -ne 1 ]; then
    echo "usage: `basename $0` [Directory]"
    exit 1
fi

DirName=$1

if [ ! -d $DirName ]; then
    echo "No such a directory '$DirName'"
    exit 1
fi

Dirs=`ls $DirName`

mkdir -p $DirName.tar/zzz.tar.log
mkdir -p $DirName.tar/check
cd $DirName

for dir in $Dirs; do    
    tar cvf ../$DirName.tar/$dir.tar $dir > ../$DirName.tar/zzz.tar.log/$dir.tar_c.log 2>&1
    
    ExitCode=$?
    if [ $ExitCode -ne 0 ]; then
	echo "@@@ TAR ERROR: $ExitCode" >> ../$DirName.tar/zzz.tar.log/$dir.tar_c.log 2>&1
    fi
done

cd ../$DirName.tar

md5sum *.tar > check/cksum.org.tmp
ls -l *.tar > check/fsize.org.tmp

sort -k 2 check/cksum.org.tmp > check/cksum.org
sort -k 9 check/fsize.org.tmp > check/fsize.org

rm check/fsize.org.tmp check/cksum.org.tmp
