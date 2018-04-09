#!/bin/sh

if [ $# -eq 0 ]; then
    echo "usage: `basename $0` [tar dir] ..."
    exit 1
fi

base=`pwd`
CompCol=$HOME/analysis_code-set/scripts/convenience/CompareColumn.py

for DirName in $@; do
    cd $base/$DirName > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        continue
    fi
    
    if [ ! -e check/cksum.org ]; then
	echo "$DirName/check/cksum.org is not exist."
	continue
    fi
    if [ ! -e check/fsize.org ]; then
	echo "$DirName/check/fsize.org is not exist."
	continue
    fi
    
    gfcksum *.tar > check/cksum.jldg.tmp
    ls -l *.tar > check/fsize.jldg.tmp
    
    sort -k 4 check/cksum.jldg.tmp > check/cksum.jldg
    sort -k 9 check/fsize.jldg.tmp > check/fsize.jldg
    
    rm check/fsize.jldg.tmp check/cksum.jldg.tmp
    
    $CompCol 5 check/fsize.jldg check/fsize.org
    if [ $? -ne 0 ]; then
	continue
    fi
    $CompCol 1 check/cksum.jldg check/cksum.org
    if [ $? -ne 0 ]; then
	continue
    fi
    
    rm check/fsize.jldg check/cksum.jldg
    echo "$DirName: correct"
done
