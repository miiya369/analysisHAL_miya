#!/bin/bash

PWD=`pwd`
if [ `basename $PWD` != "analysisHAL_miya" ]; then
    echo "This script should be carry out on the top of analysisHAL_miya/"; exit
fi

git_home=https://github.com/miiya369/analysisHAL_miya.git
git_branch=develop_miiya

git checkout $git_branch

make clean
./scripts/new/conveniences/clean.sh

current_version_all=`git log | grep Version | head -n 1`

current_version_main=`echo $current_version_all | cut -d " " -f 2 | cut -d "." -f 1`
current_version_sub=` echo $current_version_all | cut -d " " -f 2 | cut -d "." -f 2`

next_version_sub=$((current_version_sub + 1))

git add .

git commit -m "Version $current_version_main.$next_version_sub"

git push $git_home $git_branch
