#!/bin/bash

doi_san_script=/gwfefs/data/G16030/miyamoto/work.K/hadron_force/hadron_force.2016-07-06/tests/test.hadron_force.2b/scripts

if [ ! -e results ]; then
    echo "The results directory is not found."
    exit -1
fi

ln -s $doi_san_script scripts

./scripts/split.corr.csh

unlink scripts
