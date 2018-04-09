#!/bin/bash

PWD=`pwd`
if [ `basename $PWD` != "analysisHAL_miya" ]; then
    echo "This script should be carry out on the top of analysisHAL_miya/"; exit
fi

echo -n "Release the new version? (yes or other) "; read yn
if [ " $yn" != " yes" ]; then echo "Exit"; exit; fi

git_home=https://github.com/miiya369/analysisHAL_miya.git
git_branch=develop_miiya

git checkout $git_branch

cat include/Makefile.inc | sed -e "s/g++-5/g++/g" > include/Makefile.inc.tmp
mv include/Makefile.inc.tmp include/Makefile.inc

make clean
./scripts/new/conveniences/clean.sh

current_version_all=`git log | grep Version | head -n 1`

current_version_main=`echo $current_version_all | cut -d " " -f 2 | cut -d "." -f 1`

next_version_main=$((current_version_main + 1))

git add .

git commit -m "Version $current_version_main.release"

git push $git_home $git_branch

git checkout master

git merge $git_branch # We may make commit name and some comments here

git branch --delete $git_branch

git push $git_home master

git branch $git_branch

git checkout $git_branch

cat include/Makefile.inc | sed -e "s/g++/g++-5/g" > include/Makefile.inc.tmp
mv include/Makefile.inc.tmp include/Makefile.inc

git add .

git commit -m "Version $next_version_main.0"

git push $git_home $git_branch
