#!/bin/bash

tmp_dir=/bgdata4/scnfqcd/miyamoto/tmp
doi_san_script=/bghome/scnfqcd/tkymiya/hf/hadron_force.2016-01-07/tests/test.hadron_force.2b/scripts

if [ $# -ne 1 ]; then
    echo "usage: sh `basename $0` [Results Directory]"; exit -1
fi

if [ `echo $1 | cut -c 1` == '/' ]; then
    ResultsDir=$1
else
    ResultsDir=`pwd`/$1
fi

count=0
while [ -d $tmp_dir/tmp_$count ]; do
    count=$((count+1))
done
mkdir -p $tmp_dir/tmp_$count
cd       $tmp_dir/tmp_$count

ln -s $doi_san_script scripts
ln -s $ResultsDir     results

./scripts/split.corr.csh && \
mkdir results/no_split && \
mv results/correlator.* results/no_split && \
mv results.split/* results && \
rm -r results.split

unlink scripts
unlink results

cd ../
rm -rf tmp_$count
