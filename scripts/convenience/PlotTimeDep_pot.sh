#!/bin/sh

if [ $# -lt 4 ]; then
    echo "usage: `basename $0` [pot_err file] [t_min] [t_max] [r1] [r2] ..."
    echo "     : At the file name, \"REPT\" will be replaced by each time."
    exit 1
fi

args=($@)
Narg=$#

ifname=${args[0]}
t_min=${args[1]}
t_max=${args[2]}

for i in `seq 3 1 $((Narg-1))`; do
    echo "# r = ${args[$i]}"
    for t in `seq $t_min 1 $t_max`; do
	IFNAME=`echo $ifname | sed -e s/REPT/$t/g`
	echo -n "$t "
	tmp=`cat $IFNAME | grep -w "${args[$i]}"`
	if [ `echo $tmp | wc -w` -ne 8 ]; then
	    echo "There are many matching. Select another r."
	    exit 1
	fi
	echo $tmp | cut -d " " -f 5-8
    done
    if [ $i -ne $((Narg-1)) ]; then
	echo
	echo
    fi
done
