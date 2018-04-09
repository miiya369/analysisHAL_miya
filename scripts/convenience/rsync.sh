#!/bin/sh

if [ $# -le 1 ]; then
    echo "usage: `basename $0` [File(s) / Directory(s)] [dst]"
    exit 1
fi

echo "@@@ rsync start : `date`"
echo
rsync -av $@
echo
echo "@@@ rsync  end  : `date`"
