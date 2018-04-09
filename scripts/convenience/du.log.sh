#!/bin/sh

if [ -f zzz.du.log ]; then
    if [ `head -n 1 zzz.du.log` != `pwd` ]; then
	echo "This log does not correspond to current directory."
	exit 1
    fi
else
    pwd   > zzz.du.log
    echo >> zzz.du.log
fi

echo "=== `date` ===" >> zzz.du.log

du -hs * >> zzz.du.log

echo >> zzz.du.log
