#!/bin/bash

if [ $# -ne 2 ]; then
    echo "usage: `basename $0` [conf_lst.txt] [#.div]"; exit -1
fi

clst=$1
Ndiv=$2

Nconf=`cat $clst | wc -l`

if [ $((Nconf%Ndiv)) -ne 0 ]; then
    echo "ERROR: #.conf(=$Nconf) % #.div(=$Ndiv) != 0, exit."; exit -1
fi

Nconf_each=$((700/Ndiv))

for i in `seq 0 1 $((Ndiv-1))`; do 
    cat $clst | head -n $(((i+1)*Nconf_each)) | tail -n $Nconf_each #> `basename $clst`.`printf %03d $i`
    echo ""
done

