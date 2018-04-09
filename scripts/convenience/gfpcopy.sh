#!/bin/sh

if [ $# -ne 3 ]
then
    echo "usage: $0 [Copy From] [Copy To] [log-file base]"
    exit 1
fi

From=$1
To=$2
Log=$3

echo "@@@ gfpcopy start: `date`" > $Log.gfpcopy.log 2>&1

gfpcopy -p $From $To >> $Log.gfpcopy.log 2>&1

ExitCode=$?
if [ $ExitCode -eq 0 ]
then
    echo "@@@ gfpcopy Successful" >> $Log.gfpcopy.log 2>&1
else
    echo "@@@ gfpcopy ERROR: $ExitCode" >> $Log.gfpcopy.log 2>&1
fi
echo "@@@ gfpcopy  end : `date`" >> $Log.gfpcopy.log 2>&1
