#!/bin/sh

if [ $# -ne 2 ]
then
    echo "usage: sh $0 [ens] [Tshift]"
    exit 1
fi

path=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/pacs-cs_ens

unlink results_s
unlink results_c

ln -s $path$1.RHQ/s_quark/Tshift`printf %03d $2` results_s
ln -s $path$1.RHQ/c_quark/Tshift`printf %03d $2` results_c
