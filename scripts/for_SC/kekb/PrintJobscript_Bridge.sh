#!/bin/sh

if [ $# -ne 7 ]; then
    echo "usage : sh $0 [Que] [Wdir] [JobName] [Ofile] [Efile] [ParamFile] [gfile]"
    exit 1
fi

que=$1
WorkDir=$2
JobName=$3
OfileName=$4
EfileName=$5
ParamFileName=$6
gfileName=$7

echo "jobclass=qb$que"
echo "inputfile="
echo "workingdir=$WorkDir"
echo "jobname=$JobName"
echo "outputfile=$OfileName"
echo "errorfile=$EfileName"
echo "executable=$WorkDir/run.x"
echo 'environment="OMP_NUM_THREADS=8:MUSPI_NUMINJFIFOS=8:MUSPI_NUMBATIDS=10"'
echo "tasks_per_node=8"
echo "mapfile="
echo "arguments=-params $ParamFileName -f $gfileName"
