#!/bin/sh

doi_san_script=/bghome/scnfqcd/tkymiya/hf/hadron_force.2016-01-07/tests/test.hadron_force.2b/scripts

if [ ! -e "./results_s" ]
then
    echo "./results_s is not found. Maybe you wrong the directory."
fi

if [ ! -e "./results_c" ]
then
    echo "./results_c is not found. Maybe you wrong the directory."
fi

echo -n "Split Corr. OK ? : "
read yn

if [ $yn != "yes" ]
then
    echo "Exit is called. If you want to done, type 'yes'."
    exit -1
fi

ln -s $doi_san_script scripts && \
ln -s results_s results && \
./scripts/split.corr.csh && \
mkdir results/no_split && \
mv results/correlator.* results/no_split && \
mv results.split/* results && \
unlink results && \
ln -s results_c results && \
./scripts/split.corr.csh && \
mkdir results/no_split && \
mv results/correlator.* results/no_split && \
mv results.split/* results && \
unlink results

unlink scripts && rm -r results.split
