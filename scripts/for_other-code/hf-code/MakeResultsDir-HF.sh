#!/bin/sh

if [ $# -lt 3 ]; then
    echo " usage : sh `basename $0` [working dir] [conf list] [dir name 1] [dir name 2] ..."
    exit 1
fi

working_dir=$1
conf_list=$2
directories=($@)

for CONF in `cat $conf_list`; do
    for i in `seq 3 1 $#`; do
	mkdir -p "$working_dir/${directories[i-1]}/$CONF"
    done
done
